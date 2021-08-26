/**
 * In August 2021 we moved the ncov documentation from github pages
 * to docs.nextstrain.org. In order to create redirects from the former
 * to the latter, we need to create simple HTML pages at each github-pages
 * page which instruct the client to redirect. This script exists to automate
 * the creation of such HTML files.
 * See https://opensource.com/article/19/7/permanently-redirect-github-pages.
 *
 * The GitHub pages were sourced from the `/docs` directory of the `master`
 * branch. This script & the resulting HTML redirects are intended to live
 * on the `redirect-tutorial-dont-delete` branch, where gh-pages are to be sourced.
 */

const fs = require("fs");
const path = require("path");

/* ------------------------------------------------------------------------ */

const mainRTD = `https://docs.nextstrain.org/en/latest`;
const ncovRTD = `https://docs.nextstrain.org/projects/ncov/en/latest`;

const redirects = [
  // MAIN DOCS SPLASH
  [`https://nextstrain.github.io/ncov/`, `${ncovRTD}`],
  // STEPS
  [`https://nextstrain.github.io/ncov/setup`, `${ncovRTD}/analysis/setup.html`],
  [`https://nextstrain.github.io/ncov/data-prep`, `${ncovRTD}/analysis/data-prep.html`],
  [`https://nextstrain.github.io/ncov/orientation-workflow`, `${ncovRTD}/analysis/orientation-workflow.html`],
  [`https://nextstrain.github.io/ncov/orientation-files`, `${ncovRTD}/analysis/orientation-files.html`],
  [`https://nextstrain.github.io/ncov/running`, `${ncovRTD}/analysis/running.html`],
  [`https://nextstrain.github.io/ncov/customizing-analysis`, `${ncovRTD}/analysis/customizing-analysis.html`],
  [`https://nextstrain.github.io/ncov/customizing-visualization`, `${ncovRTD}/analysis/customizing-visualization.html`],
  [`https://nextstrain.github.io/ncov/sharing`, `${ncovRTD}/visualization/sharing.html`],
  [`https://nextstrain.github.io/ncov/interpretation`, `${ncovRTD}/visualization/interpretation.html`],
  [`https://nextstrain.github.io/ncov/narratives`, `${ncovRTD}/visualization/narratives.html`],
  [`https://nextstrain.github.io/ncov/multiple_inputs`, `${ncovRTD}/reference/multiple_inputs.html`],
  [`https://nextstrain.github.io/ncov/configuration`, `${ncovRTD}/reference/configuration.html`],
  [`https://nextstrain.github.io/ncov/naming_clades`, `${ncovRTD}/reference/naming_clades.html`],
  [`https://nextstrain.github.io/ncov/metadata-fields`, `${ncovRTD}/reference/metadata-fields.html`],
  [`https://nextstrain.github.io/ncov/change_log`, `${ncovRTD}/reference/change_log.html`],
  [`https://nextstrain.github.io/ncov/data_submitter_faq`, `${ncovRTD}/reference/data_submitter_faq.html`],
/** these aren't in the subproject as of the time of these redirects being implemented
 * dev_docs.md
 * glossary.md
 * translation_docs.md */
];


function generateHtml(newUrl) {
  return `
    <!DOCTYPE HTML>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <meta http-equiv="refresh" content="0;url=${newUrl}" />
            <link rel="canonical" href="${newUrl}" />
        </head>
        <body>
            <h3>
              The Nextstrain SARS-CoV-2 tutorial is now hosed on <a href="https://docs.nextstrain.org">docs.nextstrain.org</a>
            </h2>
            <p>
              This particular page has been moved to <a href="${newUrl}">${newUrl}</a>
            </p>
        </body>
    </html>
  `;
}

/**
 * There are 2 HTML files for most "normal" gh-pages URL reflecting the 2 valid URLs:
 * https://nextstrain.github.io/ncov/releases/changelog and
 * https://nextstrain.github.io/ncov/releases/changelog/index.html
 */
function githubPagesUrlToFilenames(url) {
  const subFolderStruct = url.replace("https://nextstrain.github.io/ncov", "");
  const fnames = [
    path.join(__dirname, "docs", subFolderStruct) + ".html",
    path.join(__dirname, "docs", subFolderStruct, "index.html")
  ];
  if (subFolderStruct==="/") return [fnames[1]];
  return fnames;
}


/* ------------------------------------------------------------------- */
function main() {
  redirects.forEach(([oldUrl, newUrl]) => {

    githubPagesUrlToFilenames(oldUrl).forEach((fname) => {
      if (!fs.existsSync(path.dirname(fname))) {
        fs.mkdirSync(path.dirname(fname));
      }
      fs.writeFileSync(fname, generateHtml(newUrl), {encoding: "utf8"});
    });
  });
}
main();
