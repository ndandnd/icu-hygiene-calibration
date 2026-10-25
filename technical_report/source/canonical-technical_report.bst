ENTRY { canonicalsort canonicaltechnical_reportlabel canonicaltechnical_reportreference } {} {}
FUNCTION {article} { }
FUNCTION {book} { }
FUNCTION {presort} { canonicalsort 'sort.key$ := }
FUNCTION {begin.bib} { "\begin{thebibliography}{99}" write$ newline$ }
FUNCTION {output.entry} {
  "\bibitem[" write$ canonicaltechnical_reportlabel write$ "]{" write$ cite$ write$ "}" write$ newline$
  canonicaltechnical_reportreference write$ newline$ newline$
}
FUNCTION {end.bib} { "\end{thebibliography}" write$ newline$ }
READ
ITERATE {presort}
SORT
EXECUTE {begin.bib}
ITERATE {output.entry}
EXECUTE {end.bib}
