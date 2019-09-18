"use strict";

let support_container_id = "support_cards";
let support_card_class = "gv-support";
let support_card_inner_class = "gv-support-inner";

// Focus on a specific support card.
function focusSupportCard(card_id = null) {
  let card = document.getElementById(card_id);
  for (let c of document.getElementsByClassName(support_card_class)) {
    if (c !== card) {
      c.children[0].style.height = "0%";
      c.style.height = "0%";
    }
  }
  if (card_id !== null) {
    card.children[0].style.height = "100%";
    card.style.height = "100%";
  } else {
    toggleNavigation();
  }
}

// Callback function to add the response text as support card.
function setResponseAsSupportCard(card_id) {
  return function(request) {
    // Add the request content into the content item.
    document.getElementById(card_id).children[0].innerHTML =
      request.responseText;
  };
}

// Add a support card and open it.
function tryAddSupportCard(url, card_id, gain_focus = true) {
  // Always disable the menu.
  toggleNavigation();

  // Check whether the element already exists.
  let element = document.getElementById(card_id);
  if (element === null) {
    // Create a new section element in the container that indicates that the content is loading.
    document.getElementById(support_container_id).innerHTML +=
      '<section id="' +
      card_id +
      '" class="aga-panel aga-card ' +
      support_card_class +
      '" style="height: 0%;">' +
      "<div class=" +
      support_card_inner_class +
      '><div class="aga-fill aga-center"><p>Loading content...</p></div></div>' +
      "</section>";

    // Get the item content from an external URL.
    ajaxGet(url, setResponseAsSupportCard(card_id));
  }

  if (gain_focus) focusSupportCard(card_id);
}

function initMap() {
  var mymap = L.map("map").setView([51.448, 5.49], 16);

  L.tileLayer(
    "https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw",
    {
      maxZoom: 18,
      attribution:
        'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, ' +
        '<a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, ' +
        'Imagery © <a href="https://www.mapbox.com/">Mapbox</a>',
      id: "mapbox.streets"
    }
  ).addTo(mymap);

  // L.marker([51.448, 5.49]).addTo(mymap);

  /*L.circle([51.448, 5.49], {
    color: "red",
    fillColor: "#f03",
    fillOpacity: 0.5,
    radius: 500
  }).addTo(mymap);*/

  // Make sure that the navigation menu and support cards always show on top of the map.
  let map_max_z_index = Math.max(
    getComputedStyle(document.querySelector(".leaflet-control")).zIndex,
    getComputedStyle(document.querySelector(".leaflet-top")).zIndex
  );
  let stylesheet = document.styleSheets[document.styleSheets.length - 1];
  stylesheet.insertRule(".aga-card { z-index: " + (map_max_z_index + 1) + ";}");
  stylesheet.insertRule(
    ".aga-header { z-index: " + (map_max_z_index + 2) + ";}"
  );
}
