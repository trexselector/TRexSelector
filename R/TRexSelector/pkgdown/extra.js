/* The page position and the previous/next links of the guide pages on the pkgdown site, as the
   vignettes show them. The numbered pages are the v01 .. v13 articles of the navbar menu and the
   overview is the "Get started" article, so a new page needs no change here. */
document.addEventListener("DOMContentLoaded", function () {
  var here = location.href.split("#")[0].split("?")[0];
  var links = Array.prototype.slice.call(document.querySelectorAll("nav a.dropdown-item"));
  var pages = links.filter(function (a) { return /\/articles\/v\d{2}-[^/]*\.html$/.test(a.href); })
    .map(function (a, i) { return { n: i + 1, href: a.href, title: a.textContent.trim() }; });
  var homeLink = document.querySelector('nav a.nav-link[href$="articles/TRexSelector.html"]');
  var home = homeLink ? { href: homeLink.href, title: "Overview" } : null;
  var cur = -1;
  pages.forEach(function (p, i) { if (p.href === here) cur = i; });
  var atHome = home !== null && home.href === here;
  var main = document.querySelector("main");
  if (!main || pages.length === 0 || (cur < 0 && !atHome)) return;

  /* the position above the title */
  var header = main.querySelector(".page-header h1");
  if (header && cur >= 0) {
    var where = document.createElement("p");
    where.className = "trex-where";
    where.textContent = "Page " + (cur + 1) + " of " + pages.length;
    header.parentNode.insertBefore(where, header);
  }

  /* previous and next at the end of the page; page 1 leads back to the overview,
     the overview leads to page 1 */
  var link = function (cls, href, text) {
    var a = document.createElement("a");
    a.className = cls;
    a.href = href;
    a.textContent = text;
    return a;
  };
  var pager = document.createElement("nav");
  pager.className = "trex-pager";
  pager.setAttribute("aria-label", "Guide pages");
  if (cur > 0) {
    pager.appendChild(link("prev", pages[cur - 1].href, "← " + pages[cur - 1].n + ". " + pages[cur - 1].title));
  } else if (cur === 0 && home) {
    pager.appendChild(link("prev", home.href, "← " + home.title));
  }
  if (cur >= 0) {
    var pos = document.createElement("span");
    pos.className = "where";
    pos.textContent = (cur + 1) + " / " + pages.length;
    pager.appendChild(pos);
  }
  if (cur >= 0 && cur < pages.length - 1) {
    pager.appendChild(link("next", pages[cur + 1].href, pages[cur + 1].n + ". " + pages[cur + 1].title + " →"));
  } else if (atHome) {
    pager.appendChild(link("next", pages[0].href, pages[0].n + ". " + pages[0].title + " →"));
  }
  main.appendChild(pager);
});
