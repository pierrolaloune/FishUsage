# Contributing to FishUsage

This repository holds the data and code of the paper:

> Bouchet, P., Brosse, S. & Toussaint, A. Human targeting of morphologically
> unique fishes amplifies the risk of functional erosion. *Nature
> Communications* (2026). <https://doi.org/10.1038/s41467-026-78568-9>

Its purpose is to make the published analyses reproducible. The published
results are fixed. Reports of errors, problems running the code and questions
on the data are welcome, as are corrections and improvements to the code.

## Where to contribute

Use this GitHub repository for all contributions: issues to report a problem or
ask a question, pull requests to propose a change. The repository holds the
current version of the code and documentation and keeps a public record of
each discussion.

The Zenodo deposit ([10.5281/zenodo.21873314](https://doi.org/10.5281/zenodo.21873314))
is a fixed archive of the version used in the paper. It is not updated and is
not a place for contributions.

## Reporting a problem

Before opening an issue, check that the problem has not already been reported.
Then open an issue with the template that fits:

- **Bug report**: a script fails or gives unexpected output.
- **Data question**: a question or suspected error on a dataset, species or
  variable.

Include the script concerned, the full error message, your R version and
operating system, and whether you ran `script/00_RunAll.R` or the script on its
own. A short reproducible example helps.

Errors in the source databases (FISHMORPH, FishBase, IUCN Red List, PhyloPic)
should also be reported to their maintainers, since changes made here will not
reach the original sources.

## Proposing a change

1. Fork the repository and create a branch from `main` with a short,
   descriptive name.
2. Keep each pull request to one change, and explain what it changes and why.
3. Follow the conventions of the pipeline:
   - do not rename the scripts; their numbers set the order in which
     `script/00_RunAll.R` runs them;
   - match the style of the surrounding code;
   - keep the steps flagged `[LONG]` commented out, and reload their results
     from `dataPrepared/` and `output/`;
   - do not commit large result files to `output/` or `figures/` unless the
     change requires it, and say so in the pull request.
4. Check that the scripts you changed still run from the project root and that
   `script/00_make_figures.R` still draws the figures.
5. Open a pull request with the template provided.

A change that would alter a published result will be discussed in the pull
request before it is merged, and documented in the README.

## Data and licences

The repository is distributed under the
[Creative Commons Attribution-NonCommercial 4.0 International](https://creativecommons.org/licenses/by-nc/4.0/)
licence, as archived on Zenodo. The paper itself is published under the
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International](https://creativecommons.org/licenses/by-nc-nd/4.0/)
licence.

Third-party data keep the terms of their sources: FISHMORPH, FishBase, the
IUCN Red List and PhyloPic (credits in `output/phylopic_credits.csv`). Reuse of
these data must follow the terms of each source. By contributing, you agree
that your contribution is distributed under the licence of the repository.

## Citation

If you use this code or data, please cite the paper and the Zenodo archive:

- Bouchet, P., Brosse, S. & Toussaint, A. *Nature Communications* (2026).
  <https://doi.org/10.1038/s41467-026-78568-9>
- Zenodo archive: <https://doi.org/10.5281/zenodo.21873314>

## Conduct

Keep discussions courteous, factual and focused on the work.

## Contact

Pierre Bouchet, CRBE, Université de Toulouse, France
(<pierre.bouchet@utoulouse.fr>, <pierrebdef@gmail.com>)
