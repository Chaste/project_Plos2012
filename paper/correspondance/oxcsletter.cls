
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{oxcsletter}

\def\textest{pdflatex}
\DeclareOption{xelatex}{\def\textest{xelatex}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions

\LoadClass[a4paper]{letter}
\RequirePackage{graphicx}
\RequirePackage{ifthen}

\ifthenelse{\equal{\textest}{xelatex}}
{
\RequirePackage{xltxtra}
\setmainfont[Mapping=tex-text]{Arial}
\setsansfont{FoundrySterling-Book}
}
{
\RequirePackage[T1]{fontenc}
\renewcommand{\rmdefault}{phv}
\renewcommand{\sfdefault}{foundrysterlingbook}
}

\setlength{\textheight}{200mm}
\setlength{\textwidth}{160mm}
\setlength{\topmargin}{-43mm}
\setlength{\headheight}{30mm}
\setlength{\footskip}{10mm}
\setlength{\oddsidemargin}{0mm} 
\setlength{\evensidemargin}{0mm}
\setlength{\parskip}{1em}
\setlength{\parindent}{0.0em}

\newcommand*{\refcode}[1]{\def\fromrefcode{#1}}
\newcommand*{\subject}[1]{\def\fromsubject{#1}}
\newcommand*{\position}[1]{\def\fromposition{#1}}
\newcommand*{\phone}[1]{\def\fromphone{#1}}
\newcommand*{\fax}[1]{\def\fromfax{#1}}
\newcommand*{\econtact}[1]{\def\fromecontact{#1}}

\refcode{}\subject{}
\position{}\phone{}\fax{273839}\econtact{enquiries}

\renewcommand{\ps@firstpage}{
  \setlength{\headsep}{60mm}
  \setlength{\textwidth}{166mm}
  \renewcommand{\@oddhead}{
    \parbox[t]{\textwidth}{
      \raisebox{-4mm}{{\parbox[t]{100mm}{
            {\large \textsf{Computational Biology\\Department of Computer Science\\
                Wolfson Building, Parks Road, Oxford OX1 3QD}}\\[30mm]
            \toname\\
            \toaddress
          }
        }}\hfill
      \parbox[t]{60mm}{
        \flushright \includegraphics[width=32mm]{oxlogo}\\[5mm]
        \ifthenelse{\equal{\fromrefcode}{}}
        {}{Ref: \fromrefcode\hspace*{6mm}\vspace{\baselineskip}\par}
        \@date\hspace*{6mm}
      }
    }
  }

  \renewcommand{\@oddfoot}{%
    \hspace*{-3mm}
    \parbox{\textwidth}{
      \centering
       \scriptsize 
       \textsf{General Enquiries Tel: +44 (0)1865 273838 \quad 
         \ifthenelse{\equal{\fromphone}{}}
         {}{Direct: +44 (0)1865 \fromphone}\\
         Fax: +44 (0)1865 {\fromfax} \quad
         Email: {\fromecontact}@cs.ox.ac.uk \quad
         Web: http://www.cs.ox.ac.uk
         \ifthenelse{\equal{\fromecontact}{enquiries}}
         {}{/people/\fromecontact}
       }
     }
   }
  \renewcommand{\@evenhead}{}
  \renewcommand{\@evenfoot}{}
}

\renewcommand{\ps@headings}{
  \setlength{\headsep}{30mm}
  \renewcommand{\@oddhead}{
    \parbox[t]{158mm}{
      \flushright{\includegraphics[width=20mm]{oxlogo}}
    }}

  \renewcommand{\@oddfoot}{}
  \renewcommand{\@evenhead}{\@oddhead}
  \renewcommand{\@evenfoot}{\@oddfoot}
}

\providecommand{\@evenhead}{}\providecommand{\@oddhead}{}
\providecommand{\@evenfoot}{}\providecommand{\@oddfoot}{}

\pagestyle{headings}

\renewcommand{\today}{%
  \ifcase\day\or
  1st\or 2nd\or 3rd\or 4th\or 5th\or
  6th\or 7th\or 8th\or 9th\or 10th\or
  11th\or 12th\or 13th\or 14th\or 15th\or
  16th\or 17th\or 18th\or 19th\or 20th\or
  21st\or 22nd\or 23rd\or 24th\or 25th\or
  26th\or 27th\or 28th\or 29th\or 30th\or
  31st\fi~\ifcase\month\or January\or February\or March\or
  April\or May\or June\or July\or August\or September\or October\or
  November\or December\fi \space \number\year
}

\renewcommand{\opening}[1]{\thispagestyle{firstpage}
  \setlength{\textheight}{230mm}
  #1 \par
  \ifthenelse{\equal{\fromsubject}{}}
  {}{\textbf{\fromsubject} \par}
  \nobreak
}

\renewcommand{\closing}[1]{\par\nobreak\vspace{2\parskip}
  \stopbreaks
  \noindent
  \parbox{2.5\indentedwidth}{\raggedright
    \ignorespaces #1\\[3\medskipamount]
    \fromsig\\ \position}
  \par
}
