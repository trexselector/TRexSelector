/* The running top bar with the page position on every page, the previous and next page links with
   the position between them, appended to the text block of every guide page (the contents box lays
   out the text and the page's table of contents side by side), and the dressing of the navigation
   tree: the page entries carry the short titles of the list below, their section entries are hidden
   and the bibliography entry opens the look-ups with a rule. A new page is one more entry in the
   list, in reading order. The main page is the overview: page 1 leads back to it and it leads to
   page 1. */
function trexAddNav() {
  var home = ["index.html", "Overview"];
  var pages = [
    ["getting_started.html", "Getting started"],
    ["dimensionality.html", "Dimensionality and FDR"],
    ["trex_selector.html", "The T-Rex selector"],
    ["trex_da.html", "T-Rex+DA"],
    ["trex_gvs.html", "T-Rex GVS"],
    ["screening.html", "Screen-T-Rex"],
    ["trex_spca.html", "T-Rex SPCA"],
    ["trex_tikhonov.html", "T-Rex Tikhonov"],
    ["tsolvers.html", "The T-solvers"],
    ["ml_methods.html", "The ML methods"],
    ["memory_mapping.html", "Memory mapping"],
    ["no_copies.html", "Without data copies"],
    ["parallelism.html", "Parallelism"]
  ];
  var here = location.pathname.split("/").pop();
  var atHome = here === home[0] || here === "";
  var cur = -1;
  pages.forEach(function (p, i) { if (p[0] === here) cur = i; });

  /* the top bar on every page; the brand leads to the main page */
  var bar = document.createElement("header");
  bar.className = "trex-topbar";
  bar.innerHTML = '<a class="brand" href="' + home[0] + '">TRexSelector</a>' +
                  '<span class="sub">C++ library guide</span><span class="where"></span>';
  if (cur >= 0) bar.querySelector(".where").textContent = "Page " + (cur + 1) + " of " + pages.length;
  document.body.insertBefore(bar, document.body.firstChild);

  /* the tree is built by navtree.js after this script runs and grows when a node is expanded, so
     its entries are dressed whenever it changes */
  var short = {};
  pages.forEach(function (p) { short[p[0]] = p[1]; });
  var tree = document.getElementById("nav-tree-contents");
  var dress = function () {
    Array.prototype.slice.call(tree.querySelectorAll("a[href]")).forEach(function (a) {
      var file = a.getAttribute("href").split("/").pop().split("#")[0];
      var li = a.closest("li");
      if (short[file] && li) {
        var span = a.querySelector("span");
        if (span && span.textContent !== short[file]) span.textContent = short[file];
        li.classList.add("trex-page");
      }
      if (file === "citelist.html" && li) li.classList.add("trex-lookups");
    });
  };
  if (tree) {
    dress();
    new MutationObserver(dress).observe(tree, { childList: true, subtree: true });
  }

  /* previous and next with the position between them */
  if (cur < 0 && !atHome) return;
  var contents = document.querySelector("div.contents div.textblock") ||
                 document.querySelector("div.contents");
  if (!contents) return;

  var nav = document.createElement("nav");
  nav.className = "doc-pager";
  if (cur > 0) {
    nav.innerHTML += '<a class="prev" href="' + pages[cur - 1][0] + '">&larr; ' + cur + '. ' +
                     pages[cur - 1][1] + '</a>';
  } else if (cur === 0) {
    nav.innerHTML += '<a class="prev" href="' + home[0] + '">&larr; ' + home[1] + '</a>';
  }
  if (cur >= 0) nav.innerHTML += '<span class="where">' + (cur + 1) + ' / ' + pages.length + '</span>';
  if (cur >= 0 && cur < pages.length - 1) {
    nav.innerHTML += '<a class="next" href="' + pages[cur + 1][0] + '">' + (cur + 2) + '. ' +
                     pages[cur + 1][1] + ' &rarr;</a>';
  } else if (atHome) {
    nav.innerHTML += '<a class="next" href="' + pages[0][0] + '">1. ' + pages[0][1] + ' &rarr;</a>';
  }
  contents.appendChild(nav);
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", trexAddNav);
} else {
  trexAddNav();
}
