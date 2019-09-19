/*
    Copyright 2019 Netherlands eScience Center and TU Eindhoven
    Licensed under the Apache License, version 2.0. See LICENSE for details.
*/

'use strict';

const SUPPORT_CONTAINER_ID = 'support_cards';
const SUPPORT_CARD_CLASS = 'gv-support';
const SUPPORT_CARD_INNER_CLASS = 'gv-support-inner';

// Focus on a specific support card.
function focusSupportCard(cardId = null) {
  let card = document.getElementById(cardId);
  for (let c of document.getElementsByClassName(SUPPORT_CARD_CLASS)) {
    if (c !== card) {
      c.children[0].style.height = '0%';
      c.style.height = '0%';
    }
  }
  if (cardId !== null) {
    card.children[0].style.height = '100%';
    card.style.height = '100%';
  } else {
    toggleNavigation();
  }
}

// Callback function to add the response text as support card.
function setResponseAsSupportCard(cardId) {
  return function(request) {
    // Add the request content into the content item.
    document.getElementById(cardId).children[0].innerHTML =
      request.responseText;
  };
}

// Add a support card and open it.
function tryAddSupportCard(url, cardId, gainFocus = true) {
  // Always disable the menu.
  toggleNavigation();

  // Check whether the element already exists.
  const element = document.getElementById(cardId);
  if (element === null) {
    // Create a new section element in the container that indicates that the
    // content is loading.
    // Note that the initial height is set to 0 to hide the card.
    document.getElementById(SUPPORT_CONTAINER_ID).innerHTML +=
      '<section id="' +
      cardId +
      '" class="aga-panel aga-card ' +
      SUPPORT_CARD_CLASS +
      '" style="height: 0;">' +
      '<div class=' +
      SUPPORT_CARD_INNER_CLASS +
      '><div class="aga-fill aga-center"><p>Loading content...</p></div></div>' +
      '</section>';

    // Get the item content from an external URL.
    ajaxGet(url, setResponseAsSupportCard(cardId));
  }

  if (gainFocus) focusSupportCard(cardId);
}

// Initialize the map block using Leaflet.
function initMap() {
  var map = L.map('map').fitWorld();

  // Use CartoDB wihtout labels as base layer.
  getTileLayer(tileServer.CARTODB_BASE).addTo(map);

  // Add CartoDB labels as separate pane.
  map.createPane('labels');
  let paneLabels = map.getPane('labels');
  paneLabels.style.zIndex = 650; // On top of markers but below pop-ups
  paneLabels.style.pointerEvents = 'none'; // Disable mouse events on this particular pane.
  getTileLayer(tileServer.CARTODB_LABEL, 'labels').addTo(map);

  // Add a scale control.
  L.control.scale().addTo(map);

  function onLocationFound(e) {
    // Focus on the user's location.
    const radius = e.accuracy / 2;

    L.marker(e.latlng)
      .addTo(map)
      .bindPopup('You are within ' + radius + ' meters from this point')
      .openPopup();

    L.circle(e.latlng, radius).addTo(map);
  }

  function onLocationError(e) {
    // Focus on TU/e.
    map.setView([51.448, 5.49], 16);
  }

  // Try to locate the user.
  map.on('locationfound', onLocationFound);
  map.on('locationerror', onLocationError);
  map.locate({ setView: true, maxZoom: 16 });

  // L.marker([51.448, 5.49]).addTo(map);

  /*L.circle([51.448, 5.49], {
      color: "red",
      fillColor: "#f03",
      fillOpacity: 0.5,
      radius: 500
    }).addTo(map);*/

  // L.polygon([[51.509, -0.08], [51.503, -0.06], [51.51, -0.047]]).addTo(map);

  // Make sure that the navigation and support cards always show on top of the map.
  const mapZIndexMax = Math.max(
    getComputedStyle(document.querySelector('.leaflet-control')).zIndex,
    getComputedStyle(document.querySelector('.leaflet-top')).zIndex
  );
  let stylesheet = document.styleSheets[document.styleSheets.length - 1];
  stylesheet.insertRule('.aga-card { z-index: ' + (mapZIndexMax + 1) + ';}');
  stylesheet.insertRule('.aga-header { z-index: ' + (mapZIndexMax + 2) + ';}');
}

const tileServer = {
  ARCGIS: 'ArcGIS world topo',
  CARTODB: 'CartoDB light',
  CARTODB_BASE: 'CartoDB light (no labels)',
  CARTODB_LABEL: 'CartoDB light (labels only)',
  MAPBOX: 'MapBox streets v4',
  MAPBOX_TERRAIN: 'MapBox terrain v2',
  OSM: 'OpenStreetMap',
  STAMEN: 'Stamen Design'
};

// Add a base map from the above selection of tile servers.
// If the pane parameter is not null, the layer is assigned to this pane.
function getTileLayer(server = tileServer.MAPBOX, pane = null) {
  let serverUrl = '';
  let tileOptions = { maxZoom: 18 };
  if (pane !== null) tileOptions.pane = pane;

  switch (server) {
    case tileServer.ARCGIS:
      serverUrl =
        'http://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}';
      tileOptions.attribution =
        'Map tiles © <a href="https://www.arcgis.com">ArcGIS</a> contributors. ' +
        'Map data © <a href="http://osm.org/copyright">OpenStreetMap</a> contributors, ' +
        'under <a href="http://creativecommons.org/licenses/by-sa/3.0">CC BY SA</a>. ' +
        'The GIS User Community.'; // Note that this attribution should actually be collected from the copyrightText field in http://services.arcgisonline.com/arcgis/rest/services/World_Topo_Map/MapServer?f=pjson
      break;
    case tileServer.CARTODB:
      // Alternative URL: 'https://{s}.basemaps.cartocdn.com/{id}/{z}/{x}/{y}.png'
      serverUrl =
        'https://cartodb-basemaps-{s}.global.ssl.fastly.net/{id}/{z}/{x}/{y}.png';
      tileOptions.attribution =
        'Map tiles from <a href="http://cartodb.com/attributions">CartoDB</a>, ' +
        'under <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>. ' +
        'Map data © <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, ' +
        'under <a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>.';
      tileOptions.id = 'light_all';
      break;
    case tileServer.CARTODB_BASE:
      serverUrl =
        'https://cartodb-basemaps-{s}.global.ssl.fastly.net/{id}/{z}/{x}/{y}.png';
      tileOptions.attribution =
        'Map tiles from <a href="http://cartodb.com/attributions">CartoDB</a>, ' +
        'under <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>. ' +
        'Map data © <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, ' +
        'under <a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>.';
      tileOptions.id = 'light_nolabels';
      break;
    case tileServer.CARTODB_LABEL:
      serverUrl =
        'https://cartodb-basemaps-{s}.global.ssl.fastly.net/{id}/{z}/{x}/{y}.png';
      tileOptions.attribution =
        'Map tiles from <a href="http://cartodb.com/attributions">CartoDB</a>, ' +
        'under <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>. ' +
        'Map data © <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, ' +
        'under <a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>.';
      tileOptions.id = 'light_only_labels';
      break;
    case tileServer.MAPBOX:
      serverUrl =
        'https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw';
      tileOptions.attribution =
        'Map tiles © <a href="https://www.mapbox.com/about/maps/">Mapbox</a>. ' +
        'Map data © <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, ' +
        'under <a href="http://creativecommons.org/licenses/by-sa/3.0">CC BY SA</a>. ' +
        '<strong><a href="https://www.mapbox.com/map-feedback/" target="_blank">Improve this map</a></strong>';
      tileOptions.id = 'mapbox.streets';
      break;
    case tileServer.MAPBOX_TERRAIN:
      serverUrl =
        'https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw';
      tileOptions.attribution =
        'Map tiles © <a href="https://www.mapbox.com/about/maps/">Mapbox</a>. ' +
        'Map data © <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, ' +
        'under <a href="http://creativecommons.org/licenses/by-sa/3.0">CC BY SA</a>. ' +
        '<strong><a href="https://www.mapbox.com/map-feedback/" target="_blank">Improve this map</a></strong>';
      tileOptions.id = 'mapbox.mapbox-terrain-v2';
      break;
    case tileServer.OSM:
      serverUrl = 'http://{s}.tile.osm.org/{z}/{x}/{y}.png';
      tileOptions.attribution =
        'Map data © <a href="http://osm.org/copyright">OpenStreetMap</a> contributors, ' +
        'under <a href="http://creativecommons.org/licenses/by-sa/3.0">CC BY SA</a>.';
      break;
    case tileServer.STAMEN:
      serverUrl = 'http://{s}.tile.stamen.com/{id}/{z}/{x}/{y}.jpg';
      tileOptions.attribution =
        'Map tiles © <a href="http://stamen.com">Stamen Design</a>, ' +
        'under <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>. ' +
        'Map data © <a href="http://openstreetmap.org">OpenStreetMap</a>, ' +
        'under <a href="http://creativecommons.org/licenses/by-sa/3.0">CC BY SA</a>.';
      tileOptions.id = 'watercolor';
      break;
    default:
      return null;
  }

  return L.tileLayer(serverUrl, tileOptions);
}
