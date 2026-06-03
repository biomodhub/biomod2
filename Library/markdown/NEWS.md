# CHANGES IN markdown VERSION 2.0

- The core function `mark()` is a thin wrapper of `litedown::mark()` now. Users are recommended to call **litedown** directly instead of through the wrapper.

# CHANGES IN markdown VERSION 1.13

- Cleaned `sourcepos` records when they come from metadata (thanks, @dmurdoch, #111).

- The **markdown** package is in the maintenance-only mode now. It is feature-complete, and will receive no updates except for fixing CRAN problems. New development will continue only in **litedown**: <https://github.com/yihui/litedown>.

# CHANGES IN markdown VERSION 1.12

- Provided three internal functions `html_document`, `html_vignette`, and `pdf_document` as compatibility layers to functions of the same names in the **rmarkdown** package (thanks, @jangorecki, #108).

- The default HTML template no longer wraps meta variables `include-before` and `include-after` inside `<div></div>`, because their values may contain incomplete HTML tags, e.g., `include-before = '<div>'` and `include-after = '</div>'`.

# CHANGES IN markdown VERSION 1.11

- Verbatim code blocks of the form ```` ```{lang attr1 attr2 ...} ```` were not correctly rendered.

# CHANGES IN markdown VERSION 1.10

- Raw blocks (```` ```{=lang} ````) were broken in the previous version when the support for code block attributes was added.

# CHANGES IN markdown VERSION 1.9

- Added support for attributes on fenced code blocks, e.g., ```` ```{.lang .class2 #id attr="value"}```` (thanks, @thothal, #106).

- Fixed the bug that the option `number_sections: true` doesn't work for HTML output when then input contains certain Unicode characters (thanks, @fyuniv, #104).

- Added support for rendering HTML Widgets such as **ggplotly** (thanks, @fyuniv, #105).

# CHANGES IN markdown VERSION 1.8

- Fixed the superfluous warning about path lengths in `mark_html()` (thanks, @kenjisato, #103).

# CHANGES IN markdown VERSION 1.7

- The `file` argument of `mark()` will be treated as a file path only if the file exists and the value is not wrapped in `I()`. Previously, it would be treated as a file path when it has a file extension, which could lead to confusing errors like #100 (thanks, @LukasWallrich).

- When there are emojis in the text, `mark()` may fail to identify and embed web resources (thanks, @tdhock, yihui/knitr#2254).

# CHANGES IN markdown VERSION 1.6

- Added support for footnotes, fenced `Div`s, section numbers, `{}` attributes for images/headings/fenced `Div`s, and appendices. See `vignette('intro', package = 'markdown')` for details.

- A lot of enhancements to the HTML slides format. See `vignette('slides', package = 'markdown')` for details.

- Added `vignette('article', package = 'markdown')` to demonstrate how to write an HTML article.

- If the input to `mark()` is a file, the output will also be a file by default. Previously the output would be text. If you want `mark()` to return text output when the input is a file, you may specify the argument `output = NULL`.

- The Markdown option `base64_images` has been renamed to `embed_resources`. This option can take two possible values, `"local"` and `"https"`, meaning whether to embed local and/or web (https) resources. You can specify none, either, or both of them. See `vignette('intro', package = 'markdown')` for details.

- Removed the option `standalone` from the list of Markdown options. Please use the argument `template = TRUE/FALSE` of `mark()` instead. The option `standalone = TRUE` was equivalent to `template = TRUE`.

- Added the option `auto_identifiers` (enabled by default) to automatically add IDs to headings, e.g., `# Hello world!` will be converted to `<h1 id="hello-world">Hello world!</h1>`. You can certainly override the automatic ID by providing an ID manually via the `{#id}` attribute, e.g., `# Hello world! {#hello}`.

- Renamed the `mathjax` option to `js_math` to allow for other JS math libraries. The default library was changed from MathJax to KaTeX. To continue using MathJax, you may set `js_math: mathjax`.

- Removed the option `mathjax_embed` from the list of Markdown options. To embed the MathJax library, enable `"https"` in the `embed_resources` option instead. Note that only MathJax v3 can be partially embedded, and lower versions cannot.

- Renamed the option `highlight_code` to `js_highlight`, and added support for an alternative syntax highlighting JS library Prism.js, which became the default. To continue using the old default `highlight.js`, you may set the `js_highlight` option to `highlight`.

- The default version of MathJax has been changed from v2 to v3.

- The default version of highlight.js has been changed from 11.6.0 to 11.7.0, and the default style has been switched from `github` to `xcode`.

# CHANGES IN markdown VERSION 1.5

- Values of meta variables `title`, `author`, and `date` (if provided) will be transformed to the target output format before they are passed into templates.

- Fixed the bug that the default CSS was not added to HTML output.

- Removed dependency on the **mime** package.

- Added experimental support for HTML slides: `markdown::mark_html(..., meta = list(css = c('default', 'slides'), js = 'slides'))`. If you prefer knitting `Rmd` documents in RStudio, you may use the output format:

  ```yaml
  output:
    markdown::html_format:
      meta:
        css: [default, slides]
        js: [slides]
  ```

  See https://yihui.org/en/2023/01/minimal-r-markdown/ for a demo.

# CHANGES IN markdown VERSION 1.4

- Empty `\title{}` in LaTeX output will be removed (along with `\maketitle`).

- highlight.js is loaded from https://www.jsdelivr.com/package/gh/highlightjs/cdn-release by default now. This means more languages are supported (not only R), but also means syntax-highlighting will not work offline at the moment (it will be improved in future).

- MathJax failed to load in the previous version. The bug has been fixed now.

- Removed the function `markdownExtensions()`.

# CHANGES IN markdown VERSION 1.3

- Switched the underlying Markdown rendering engine from the C library **sundown** (which has been deprecated for a decade) to the R package **commonmark** (thanks, @jeroen, yihui/knitr#1329).

- The functions `renderMarkdown()` and `markdownToHTML()` have been renamed to `mark()` and `mark_html()`, respectively. The old names are still kept in this package for backward-compatibility.

- Removed the arguments `stylesheet` and `fragment.only` in `mark_html()`. For `stylesheet`, please use the argument `meta = list(css = ...)` to provide the CSS stylesheet. For `fragment.only`, please use `mark_html(template = FALSE)` or `mark_html(options = '-standalone')` instead of `fragment.only = TRUE`. Currently these old arguments are still accepted internally, but may be deprecated and dropped in the long run.

- The `file` argument of `mark()` and `mark_html()` can also take a character vector of Markdown text now.

- Removed functions `rendererExists()`, `rendererOutputType()`, and `registeredRenderer()`. They were primarily for internal use.

- Deprecated the function `markdownExtensions()`. All extensions should be specified via the `options` argument of functions like `mark()`, e.g., `mark(options = '+table+tasklist')`. See all options on the help page `?markdown::markdown_options`.

- Renamed `markdownHTMLOptions()` to `markdown_options()`.

# CHANGES IN markdown VERSION 1.2

- Fixed the warnings "a function declaration without a prototype is deprecated in all versions of C" (#94).

# CHANGES IN markdown VERSION 1.1

## MAJOR CHANGES

- renderMarkdown() and markdownToHTML() will signal an error if the input file is not encoded in "UTF-8".

# CHANGES IN markdown VERSION 1.0

## MAJOR CHANGES

- The default value of the encoding argument of renderMarkdown() and markdownToHTML() has been changed from getOption("encoding") to "UTF-8". The encoding of the input file will always be assumed to be UTF-8.

- markdownToHTML() will return a character vector encoded in UTF-8 (instead of the system's native encoding) when not writing to an output file.

# CHANGES IN markdown VERSION 0.9

## BUG FIXES

- Fixed clang-UBSAN and valgrind issues (thanks, @yixuan, #92).

# CHANGES IN markdown VERSION 0.8

## MINOR CHANGES

- the MathJax CDN URL was replaced by https://www.bootcdn.cn/mathjax/

## BUG FIXES

- fixed https://github.com/rstudio/htmltools/issues/30: markdownToHTML() did not work with empty files (thanks, @VermillionAzure)

# CHANGES IN markdown VERSION 0.7.7

## BUG FIXES

- renderMarkdown() works now even if text = character(0) or ""

- added an `encoding` argument to renderMarkdown() since multi-byte characters in renderMarkdown() did not work on Windows (thanks, Kohske Takahashi, #63)

- fixed #64: invalid 'n' argument in rpubsUpload() (thanks, Wouter van Atteveldt)

## MAJOR CHANGES

- if renderMarkdown() returns a character vector, it will be marked with the UTF-8 encoding if it contains multi-byte characters

# CHANGES IN markdown VERSION 0.7.4

## NEW FEATURES

- when an image is the only element of its parent node in the HTML output document, it is automatically centered on the page

## MINOR CHANGES

- images that have already been base64 encoded will not be encoded again (#61)

- the URL of the MathJax CDN was updated to cdn.mathjax.org

# CHANGES IN markdown VERSION 0.7.2

## BUG FIXES

- fixed #60: MathJax may be included even if it is unnecessary when syntax highlighting is enabled (thanks, @aoles)

- fixed a bug which may hang R when building R Markdown vignettes in a wrong locale (thanks, Dan Tenenbaum, yihui/knitr#782)

# CHANGES IN markdown VERSION 0.7

## BUG FIXES

- if both the 'file' and 'text' arguments are provided but file = NULL, e.g. markdownToHTML(file = NULL, text = ?), markdownToHTML() can throw an error indicating the file is invalid (thanks, Tyler Rinker, hadley/staticdocs#66)

- markdownToHTML(text = ?, output = ?) was broken (#54)

# CHANGES IN markdown VERSION 0.6.5

## NEW FEATURES

- added an argument 'encoding' to markdownToHTML() to specify the character encoding of the input markdown file, and the HTML output file is always encoded in UTF-8 now (thanks, Kohske Takahashi, #50)

# CHANGES IN markdown VERSION 0.6.4

## NEW FEATURES

- added 'mathjax_embed' to HTML options for embedding the MathJax JavaScript in the HTML document rather than linking to it online. Note the JavaScript code is read from the http instead of https MathJax URL. Contributed by Henrik Bengtsson.

- added another vignette to show the HTML output of the original vignette (see browseVignettes('markdown'))

- the default CSS style was tweaked (major changes include: page width is at most 800px, more line height, slightly larger fonts, and a different syntax highlighting theme)

# CHANGES IN markdown VERSION 0.6.3

## NEW FEATURES

- added a new argument 'template' to markdownToHTML() so that we can customize the HTML template (by default, it uses the template 'resources/markdown.html' in this package); thanks, Nacho Caballero

- the options markdown.HTML.stylesheet and markdown.HTML.header used in markdownToHTML() can be character vectors (they will be processed by paste(x, collapse = '\n')

## MAJOR CHANGES

- the 'text' argument in markdownToHTML() and renderMarkdown() is treated as lines of input now, i.e. if 'text' is provided, it is passed to the markdown renderer as paste(text, collapse = '\n'); in the previous versions, it was processed by paste(text, collapse = '')

# CHANGES IN markdown VERSION 0.6

## DOCUMENTATION

- added a package vignette; see browseVignettes(package = 'markdown')

# CHANGES IN markdown VERSION 0.5.5

## NEW FEATURES

- added a new argument 'header' to markdownToHTML() to insert code into the HTML header (e.g. custom CSS styles)

## BUG FIXES

- fixed #25 and #27: minor documentation problems

- fixed #26: the HTML output file will be written relative to the current working directory now when it contains images that need to be base64 encoded

- fixed #28: the image URL should be decoded before the image is based64 encoded

## MISC

- Yihui Xie has taken over the maintainership for this package from Jeffrey Horner

# CHANGES IN markdown VERSION 0.5.4

## NEW FEATURES

- Both Pandoc title blocks and Jekyll front matter sections are skipped when rendering markdown documents.

# CHANGES IN markdown VERSION 0.5.3

## NEW FEATURES

- C/C++ is now a supported language for code block highlighting.

## MAJOR CHANGES

- 'hard_wrap' has been dropped while 'mathjax' and 'highlight_code' have been added to the default list of html options.

## BUG FIXES

- fixed parsing of math equations when located at the end of a line.

# CHANGES IN markdown VERSION 0.5.2

## NEW FEATURES

- with the new 'latex_math' markdown extensions, users can include math equations using several syntaxes. For block level equations, use $$latex ... $$, $$ ... $$, or \[ ... \]. For inline equations, use $latex...$, $...$, or \( ... \).

## MAJOR CHANGES

- the markdown extension 'ignore_math' was replaced with 'latex_math'.

- users can now use the markdown.HTML.stylesheet option to override the package default stylesheet.

- setting the fragment_only rendering option or the fragment.only parameter to markdownToHTML will base64 encode images if applicable. version 0.5.1 did not.

# CHANGES IN markdown VERSION 0.5.1

## BUG FIXES

- fixed a GUIDgenerator bug; for escaping math equations before markdown parsing begins.

- image encoding was fixed for the case when there are more than one included in a markdown document.

# CHANGES IN markdown VERSION 0.5

## NEW FEATURES

- added fragment.only parameter to markdownToHTML

- added new html rendering options base64_images, fragment_only, mathjax, and highlight_code

- added new markdown extension ignore_math

## MAJOR CHANGES

- removed safelink from default html rendering options

- the default html rendering options are now hard_wrap, use_xhtml, smartypants, and base64_images.

## BUG FIXES

- fixed syntax errors in C exports

# CHANGES IN markdown VERSION 0.4

## NEW FEATURES

- added support for post-processing HTML using smartypants filter

- added optional support for rendering a table of contents

## MAJOR CHANGES

- changed exported C functions to use an rmd_ prefix (eliminating potential symbol conflicts with other packages)

- changed default html rendering options to use_xhtml, hard_wrap, safelink, and smartypants

## BUG FIXES

- eliminated name collision with render_markdown function in knitr

