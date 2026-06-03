\name{NEWS}
\title{NEWS file for the selectr package}
\encoding{UTF-8}

\section{Changes in version 0.5-1}{
  \subsection{BUG FIXES}{
    \itemize{

      \item Fixed handling of CSS unicode escapes in attribute values.
      This would be observed when the attribute value contained hexadecimal
      sequences like \code{(abcdef)} where the characters inside the parentheses
      were not properly escaped. This fix ensures that such sequences are correctly
      translated to their XPath equivalents. Thanks to André Veríssimo for
      reporting the issue.

    }
  }
}

\section{Changes in version 0.5-0}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item Added support for CSS Selectors Level 4 pseudo-classes
      \code{:is()}, \code{:where()}, and \code{:has()}. The \code{:is()}
      pseudo-class matches elements against a list of selectors, taking the
      maximum specificity from its arguments. The \code{:where()} pseudo-class
      works similarly to \code{:is()} but always has zero specificity. The
      \code{:has()} pseudo-class represents an element if any of the relative
      selectors match when anchored against that element, with specificity
      calculated from the maximum of its arguments.

      \item Added support for complex selectors in \code{:nth-child()} and
      \code{:nth-last-child()} using the \code{of S} syntax
      (e.g., \code{:nth-child(2 of .foo)}). This allows matching the nth
      child that matches a specific selector or selector list.

      \item Extended \code{:not()} to accept multiple selectors separated by
      commas (e.g., \code{:not(.foo, #bar)}), following CSS Selectors Level 4.
      Specificity is now calculated as the maximum specificity among all
      arguments, rather than the sum.

      \item Added support for additional CSS Selectors Level 4 pseudo-classes:
      \code{:any-link}, \code{:target-within}, and \code{:local-link}. These
      pseudo-classes do not match any elements in static XML/HTML documents
      and translate to XPath expressions that always evaluate to false.
      For now, most of the new Level 4 pseudo-classes that depend on dynamic
      document state e.g. \code{:user-valid} and \code{:placeholder-shown} are
      not implemented, but may be at a future date to be non-matching selectors.

      \item The \code{:lang()} and \code{:dir()} pseudo-classes now support
      multiple comma-separated arguments (e.g., \code{:lang(en, fr, de)}).

      \item Added \code{:matches()} as a backwards-compatible alias for
      \code{:is()}.

    }
  }
  \subsection{MINOR CHANGES}{
    \itemize{

      \item Improved sibling selector translation to use a more compact form.
      For the adjacent sibling combinator \code{a + b}, the generated XPath now
      uses \code{a/following-sibling::*[1][self::b]} instead of
      \code{a/following-sibling::*[(name() = 'b') and (position() = 1)]}.

      \item The descendant combinator \code{a b} now uses
      \code{a//b} instead of \code{a/descendant::b} for a more concise XPath.
      Unfortunately a similar optimisation cannot be applied to in general
      when attempting to replace \code{descendant-or-self::} with \code{.//a}
      as it would prevent root nodes being matched correctly.

      \item Improved validation of CSS selector arguments. Better error
      messages are now provided when pseudo-elements appear inside functional
      pseudo-classes where they are not permitted (e.g., inside \code{:is()},
      \code{:matches()}, \code{:where()}, or \code{:has()}).

      \item Enhanced input validation for \code{:lang()} and \code{:dir()}
      pseudo-classes to ensure proper argument formatting and to reject
      invalid or empty language tags.

      \item Improved handling of edge cases in selector parsing, including
      better validation of class selector syntax and more robust handling of
      null or missing element components.

      \item Simplified method registration for XML and xml2 objects. No longer
      necessary to hook into package load/unload events.

    }
  }
}

\section{Changes in version 0.4-2}{
  \subsection{MINOR CHANGES}{
    \itemize{

      \item Improve handling of vectors of length > 1 in logical comparison.
      Contributed by Garrick Aden-Buie.

      \item Minor improvements to error message construction.
      Contributed by Michael Chirico.

    }
  }
}

\section{Changes in version 0.4-1}{
  \subsection{BUG FIXES}{
    \itemize{

      \item When the \pkg{R.oo} package is attached, the use of class
      selectors no longer worked. This is due to the use of the \code{Class}
      name for \pkg{R.oo}'s base class object, where \pkg{selectr} was also
      using it (but not exporting) the same name of \code{Class} for
      representing a class selector. Consequently, \pkg{selectr}'s code was
      changed to rename the class to avoid any clashing. Because it was not
      exported, this is purely an internal change. Thanks to
      Francois Lemaire-Sicre for reporting the issue.

    }
  }
}

\section{Changes in version 0.4-0}{
  \subsection{MINOR CHANGES}{
     \itemize{

      \item Large rewrite of internals to use the R6 OO system instead of
      Reference Classes. This does not affect any external facing code as the
      results should be identical to the previous implementation, which is why
      this change is marked as minor. Initial and crude performance testing
      (by running the test suite) indicates that the R6 implementation is
      approximately twice as fast at generating XPath as the Reference Classes
      implementation.

      \item The minimum required version of R for \pkg{selectr} has been
      increased from \code{2.15.2} to \code{3.0} as that is the minimum
      required version of \pkg{R6}.

      \item Minor performance enhancements have been made. Not only is
      \pkg{R6} faster than Reference Classes, the use of string formatting
      has been replaced with string concatenation. Additionally dynamic
      calling of methods via \code{do.call()} has been replaced with direct
      method calls.

     }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item The issues in previous releases where methods can sometimes be
      missing should now be resolved. The bug appeared to lie in core
      Reference Classes code. By switching to \pkg{R6}, this type of issue
      should no longer be possible.

    }
  }
}


\section{Changes in version 0.3-2}{
  \subsection{MINOR CHANGES}{
     \itemize{

      \item Improved method registration for \pkg{XML} and \pkg{xml2} 
      objects. Avoids checks on each use and is only performed once per
      dependent package load/unload.

     }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item In some environments, reference class methods were missing at
      runtime. This appears to be due to some internal behaviour in them
      \pkg{methods} package where methods are registered on an objects when
      the \code{$} operator is used for a field or method. Instead, when
      a method is missing, they are manually bound to the object.

    }
  }
}

\section{Changes in version 0.3-1}{
  \subsection{MINOR CHANGES}{
     \itemize{

      \item Enabled partial matching on the translator argument to
      \code{css_to_xpath()}. Instead of defaulting to a generic translator,
      a non-matching argument will be returned with an error.

      \item Introduced many more unit tests via the \pkg{covr} package.
      This enabled dead code to be trimmed and also identified areas of code
      which needed improvement. Minor enhancements include: tolerate
      whitespace within a \code{:not()}, more consistent results returned
      from parser methods, improvements to argument parsing.

     }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item The \code{|=} attribute matching operator was not being parsed
      correctly for the generic translator.

      \item Handle scenario where a CSS comment is unclosed. Results in
      everything after the comment start to be removed (which may or may
      not result in a valid selector).

    }
  }
}

\section{Changes in version 0.3-0}{
  \subsection{MAJOR CHANGES}{
    \itemize{

      \item Added support for documents from the \pkg{xml2} package.

      \item selectr now also does not strictly depend on the XML
      package. If either the \pkg{XML} or \pkg{xml2} packages are
      present (which are required for the \code{querySelector} methods
      to work) then \code{querySelector} will begin to work for them.
      This also enables selectr to be used for translation-only.

    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item Improve support for nth-*(an+b) selectors. Ported from
      cssselect.

    }
  }
}

\section{Changes in version 0.2-3}{
  \subsection{MINOR CHANGES}{
    \itemize{

      \item Code cleanup contributed by Kun Ren (#1).

      \item Updated DESCRIPTION to include URL and BugReports fields.
      Also update email address.

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item Fix behaviour for nth-*(an+b) pseudo-class selectors for
      negative a's. Contributed to cssselect by Paul Tremberth,
      ported to R.

      \item Escape delimiting characters to support new version of the
      stringr package. Probably should have been done in the first place.
      Reported by Hadley Wickham (#5).

    }
  }
}

\section{Changes in version 0.2-2}{
  \subsection{MINOR CHANGES}{
    \itemize{

      \item Corrected licence to BSD 3 clause. This was the licence in
      use previously, but has now been made more explicit.

      \item Removed 'Enhances' field because we import functions from
      \pkg{XML}. This choice is made because \pkg{XML} is a required
      package, rather than an optional package that can be worked with.
      This and the previous change have been made to keep up with recent
      changes in R-devel.

    }
  }
}

\section{Changes in version 0.2-1}{
  \subsection{MINOR FEATURES}{
    \itemize{

      \item Added a 'CITATION' file which cites a technical report on
      the package.

      \item \code{show()} methods are now available on internal objects,
      making interactive extensibility and bug-fixing easier. This is
      simply wrapping the \code{repr()} methods (mirroring the Python
      source) that the same objects have.

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item Use the session character encoding to determine whether to
      run unicode tests. Tests break in non-unicode sessions otherwise.

    }
  }
}

\section{Changes in version 0.2-0}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item Introduced new functions \code{querySelectorNS()} and
      \code{querySelectorAllNS()} to ease the use of namespaces within a
      document. Previously this would have required knowledge of XPath.

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item Fix meaning of \code{:empty}, whitespace is not empty.

      \item Use \code{lang()} for XML documents with the \code{:lang()}
      CSS selector.

      \item \code{|ident} no longer produces a parsing error, but is now
      equivalent to just 'ident'.

    }
  }
}

\section{Changes in version 0.1-1}{
  \subsection{BUG FIXES}{
    \itemize{

      \item Now testing unicode only in non-Windows platforms on package
      check. Output should still be consistent, just depends on the
      current charset being unicode.

    }
  }
}

\section{Changes in version 0.1-0}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item Initial port of the Python 'cssselect' package. Code is very
      literally ported, including the test suite.

      \item Wrapped translation functionality into a single function,
      \code{css_to_xpath()}.

      \item Created two convenience functions, \code{querySelector()} and
      \code{querySelectorAll()}. These mirror the behaviour of the same
      functions present in a web browser. \code{querySelector()} returns a
      node, while \code{querySelectorAll()} returns a list of nodes.

    }
  }
}
