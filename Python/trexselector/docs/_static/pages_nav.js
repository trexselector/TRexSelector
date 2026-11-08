/* The running top bar with the page position on every page, and the position "n / N" between the
   previous and next links of a guide page, as the R vignettes and the C++ guide show them. The numbered
   pages are the entries of the "User guide" list of the sidebar, so a new page needs no change here. */
document.addEventListener("DOMContentLoaded", function () {
  var list = null;
  document.querySelectorAll(".sidebar-tree > p.caption").forEach(function (p) {
    var next = p.nextElementSibling;
    if (p.textContent.trim() === "User guide" && next && next.tagName === "UL") list = next;
  });
  var items = list ? Array.prototype.slice.call(list.children).filter(function (li) {
    return li.classList.contains("toctree-l1");
  }) : [];
  var cur = -1;
  items.forEach(function (li, i) { if (li.classList.contains("current-page")) cur = i; });

  /* the top bar on every page; the brand leads to the landing page */
  var brand = document.querySelector(".sidebar-brand");
  var bar = document.createElement("header");
  bar.className = "trex-topbar";
  bar.innerHTML = '<a class="brand" href="' + (brand ? brand.getAttribute("href") : "index.html") +
                  '">trexselector</a><span class="sub">Python package guide</span>' +
                  '<span class="where"></span>';
  if (cur >= 0) bar.querySelector(".where").textContent = "Page " + (cur + 1) + " of " + items.length;
  document.body.insertBefore(bar, document.body.firstChild);

  /* the position between the previous and next links */
  var pager = document.querySelector(".related-pages");
  if (cur < 0 || !pager) return;
  var pos = document.createElement("div");
  pos.className = "page-position";
  pos.textContent = (cur + 1) + " / " + items.length;
  pager.appendChild(pos);
});
