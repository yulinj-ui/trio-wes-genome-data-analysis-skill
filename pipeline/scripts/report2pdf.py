#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report2pdf.py —— 把分析报告 .md 转成可直接邮件外发的 PDF（零外部依赖，只用本机已有组件）

用法:
    python3 report2pdf.py 报告.md                      # 输出同名 .pdf
    python3 report2pdf.py 报告.md -o /path/out.pdf
    python3 report2pdf.py 报告.md --password 8位密码    # 加密（外发患者报告建议加）
    python3 report2pdf.py 报告.md --title "左页眉" --right "右页眉"
    python3 report2pdf.py 报告.md --no-landscape       # 宽表不排横向页

技术路线（全部是本机已装的东西，不需要 brew / pandoc / LaTeX）:
    python-markdown  →  内置 print CSS  →  无头 Chrome 打印  →  PyMuPDF 加页眉页脚/加密

设计要点（都是 2026-08-26 实测踩出来的）:
  · 字体栈按"能解析就用"顺序排：PingFang SC → Heiti SC → Songti SC。
    ⚠️ 某些 Mac 未装 PingFang（本机即是），Heiti SC 可解析出 STHeitiSC-Light + Medium，
      粗体是真字重而非合成加粗；Hiragino Sans GB 在 Chrome 里解析不出来，不要放进栈里当依靠。
  · 列数 ≥8 的宽表自动单独排成 A4 横向页（否则每列被挤到 2-3 个字宽）。
  · 页眉页脚用 PyMuPDF 事后盖，因为 Chrome 无头模式不支持 CSS @page 的 @top/@bottom。
    ⚠️ 盖字必须用 fitz.Font + TextWriter 并调 subset_fonts()，
      否则每页嵌一份完整 .ttc，13 页能撑到 29MB（实测）。
"""
import sys, os, re, argparse, subprocess, tempfile, shutil

CHROME_CANDIDATES = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
]
STAMP_FONT_CANDIDATES = [
    "/System/Library/Fonts/PingFang.ttc",
    "/System/Library/Fonts/STHeiti Light.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Songti.ttc",
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
]
CJK_STACK = '"PingFang SC","Heiti SC","Songti SC","Hiragino Sans GB"'

CSS = """
@page      { size: A4;           margin: 17mm 13mm 19mm 13mm; }
@page rot  { size: A4 landscape; margin: 16mm 14mm 18mm 14mm; }
.rot { page: rot; break-before: page; break-after: page; }

html { font-size:10pt; -webkit-print-color-adjust:exact; print-color-adjust:exact; }
body { font-family:__CJK__,"Helvetica Neue",sans-serif; line-height:1.62; color:#1b1b1b; margin:0; }
h1 { font-size:19pt;   font-weight:600; color:#1d3550; margin:0 0 4mm; padding-bottom:3mm;
     border-bottom:2.2pt solid #2f4f6f; }
h2 { font-size:13.5pt; font-weight:600; color:#1d3550; margin:9mm 0 3mm; padding-bottom:1.6mm;
     border-bottom:.8pt solid #c9d4e0; break-after:avoid; page-break-after:avoid; }
h3 { font-size:11.2pt; font-weight:600; color:#2f4f6f; margin:6mm 0 2mm; break-after:avoid; page-break-after:avoid; }
h4 { font-size:10.2pt; font-weight:600; margin:4mm 0 1.5mm; break-after:avoid; }
p  { margin:0 0 2.6mm; orphans:2; widows:2; }
strong { font-weight:600; color:#111; }
hr { border:0; border-top:.6pt solid #dfe4ea; margin:6mm 0; }
ul,ol { margin:0 0 2.6mm; padding-left:6mm; }
li { margin-bottom:1.1mm; }
blockquote { margin:3mm 0; padding:2.4mm 4mm; background:#f5f8fb; border-left:2.6pt solid #5b7fa6;
             color:#22303d; break-inside:avoid; page-break-inside:avoid; }
blockquote p:last-child { margin-bottom:0; }
code { font-family:"SF Mono",Menlo,monospace; font-size:8.6pt; background:#f2f3f5;
       padding:.4mm 1mm; border-radius:2px; }
pre  { background:#f7f8fa; border:.5pt solid #e2e5e9; border-radius:3px; padding:2.6mm 3mm;
       font-size:7.4pt; line-height:1.42; overflow-wrap:anywhere; white-space:pre-wrap;
       break-inside:avoid; page-break-inside:avoid; }
pre code { background:none; padding:0; font-size:7.4pt; }
table { width:100%; border-collapse:collapse; margin:3mm 0 4mm; table-layout:fixed; }
table.norm { font-size:8.4pt; line-height:1.42; }
table.mid  { font-size:7.5pt; line-height:1.36; }
table.wide { font-size:7.6pt; line-height:1.38; }
th,td { border:.5pt solid #cfd6de; padding:1.3mm 1.6mm; vertical-align:top;
        overflow-wrap:anywhere; word-break:break-word; text-align:left; }
th { background:#eef2f7; font-weight:600; color:#1d3550; }
tr { break-inside:avoid; page-break-inside:avoid; }
thead { display:table-header-group; }
tbody tr:nth-child(even) { background:#fafbfc; }
"""


def die(msg, code=1):
    print(f"❌ {msg}", file=sys.stderr); sys.exit(code)


def find_first(paths, what):
    for p in paths:
        if os.path.exists(p):
            return p
    die(f"找不到{what}，试过：\n   " + "\n   ".join(paths))


def lint(md_text):
    """交付前排版自检：报告里最常见的两类"看起来不美观"的成因。"""
    L = md_text.split("\n"); warn = []
    for i in range(1, len(L)):
        if re.match(r"^\s*[-*+] ", L[i]) and L[i-1].strip() \
           and not re.match(r"^\s*[-*+] |^\s*\d+\. |^\s*>|^\|", L[i-1]):
            warn.append(f"第 {i+1} 行：列表前缺空行 → 会连成普通段落并原样显示 '-'")
    for i, l in enumerate(L):
        if re.match(r"^\|[\s:|-]*-[\s:|-]*\|", l):
            n = l.count("|") - 1
            if n >= 8:
                warn.append(f"第 {i} 行：表格 {n} 列（已自动排横向页；>10 列建议拆表）")
    return warn


def md_to_html(md_text, title, landscape=True):
    try:
        import markdown
    except ImportError:
        die("缺 python-markdown：pip3 install markdown")
    body = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "sane_lists", "attr_list"],
        output_format="html5",
    )

    def tag(m):
        t = m.group(0)
        ncol = t.count("<th") or max(
            (r.count("<td") for r in re.findall(r"<tr>.*?</tr>", t, re.S)), default=0)
        cls = "wide" if ncol >= 8 else ("mid" if ncol >= 6 else "norm")
        t = t.replace("<table>", f'<table class="{cls}">', 1)
        return f'<div class="rot">{t}</div>' if (cls == "wide" and landscape) else t

    body = re.sub(r"<table>.*?</table>", tag, body, flags=re.S)
    css = CSS.replace("__CJK__", CJK_STACK)
    return (f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">'
            f"<title>{title}</title><style>{css}</style></head><body>{body}</body></html>")


def stamp(raw_pdf, out_pdf, hdr_left, hdr_right, password=None):
    try:
        import fitz
    except ImportError:
        die("缺 PyMuPDF：pip3 install pymupdf")
    fitz.TOOLS.mupdf_display_errors(False)   # 屏蔽 subset_fonts 对 Chrome 字体结构的无害告警
    fp = find_first(STAMP_FONT_CANDIDATES, "中文盖字字体")
    doc = fitz.open(raw_pdf); n = len(doc); font = fitz.Font(fontfile=fp)
    for i, page in enumerate(doc):
        W, H = page.rect.width, page.rect.height
        tw = fitz.TextWriter(page.rect, color=(.54, .54, .54))
        if i:                                    # 首页不加页眉
            if hdr_left:
                tw.append((36, 30), hdr_left, font=font, fontsize=7.5)
            if hdr_right:
                tw.append((W - 36 - font.text_length(hdr_right, 7.5), 30),
                          hdr_right, font=font, fontsize=7.5)
            page.draw_line(fitz.Point(36, 34), fitz.Point(W - 36, 34),
                           color=(.85, .87, .89), width=.4)
        f = f"第 {i+1} 页 / 共 {n} 页"
        tw.append(((W - font.text_length(f, 8)) / 2, H - 26), f, font=font, fontsize=8)
        tw.write_text(page)
    doc.subset_fonts(verbose=False)              # 不做这步 13 页能撑到 29MB
    kw = dict(garbage=4, deflate=True, clean=True)
    if password:
        kw.update(encryption=fitz.PDF_ENCRYPT_AES_256,
                  user_pw=password, owner_pw=password + "_own",
                  permissions=int(fitz.PDF_PERM_PRINT | fitz.PDF_PERM_COPY |
                                  fitz.PDF_PERM_ACCESSIBILITY))
    doc.save(out_pdf, **kw); doc.close()
    return n


def main():
    ap = argparse.ArgumentParser(description="报告 .md → 可外发 PDF")
    ap.add_argument("md")
    ap.add_argument("-o", "--out")
    ap.add_argument("--title", help="左页眉；默认取 .md 的一级标题")
    ap.add_argument("--right", help="右页眉；默认自动抓 '送检号 xxx'")
    ap.add_argument("--password", help="打开密码（外发患者报告建议加，AES-256）")
    ap.add_argument("--no-landscape", action="store_true", help="宽表不排横向页")
    a = ap.parse_args()

    if not os.path.isfile(a.md):
        die(f"文件不存在: {a.md}")
    md = open(a.md, encoding="utf-8").read()
    out = a.out or os.path.splitext(a.md)[0] + ".pdf"

    h1 = re.search(r"^#\s+(.+)$", md, re.M)
    title = a.title or (h1.group(1).strip() if h1 else os.path.basename(a.md))
    m = re.search(r"(送检号|样本编号|病例号)\s*[:：]?\s*([A-Za-z0-9\-]+)", md)
    right = a.right if a.right is not None else (f"{m.group(1)} {m.group(2)}" if m else "")

    for w in lint(md):
        print(f"  ⚠️ 排版自检: {w}")

    chrome = find_first(CHROME_CANDIDATES, "Chrome/Chromium（用于打印 PDF）")
    tmp = tempfile.mkdtemp(prefix="rpt2pdf_")
    try:
        html, raw = os.path.join(tmp, "r.html"), os.path.join(tmp, "r.pdf")
        open(html, "w", encoding="utf-8").write(
            md_to_html(md, title, landscape=not a.no_landscape))
        r = subprocess.run(
            [chrome, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
             f"--print-to-pdf={raw}", "--virtual-time-budget=20000", f"file://{html}"],
            capture_output=True, text=True, timeout=300)
        if not os.path.exists(raw):
            die(f"Chrome 打印失败:\n{r.stderr[-800:]}")
        n = stamp(raw, out, title, right, a.password)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    mb = os.path.getsize(out) / 1048576
    print(f"✅ {out}\n   {n} 页 · {mb:.2f} MB"
          + ("  · 已 AES-256 加密" if a.password else ""))


if __name__ == "__main__":
    main()
