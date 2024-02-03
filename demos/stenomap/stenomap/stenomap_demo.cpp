/*
The Necklace Map console application implements the algorithmic
geo-visualization method by the same name, developed by
Bettina Speckmann and Kevin Verbeek at TU Eindhoven
(DOI: 10.1109/TVCG.2010.180 & 10.1142/S021819591550003X).
Copyright (C) 2021  Netherlands eScience Center and TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "stenomap_demo.h"

#include <QApplication>
#include <QCheckBox>
#include <QMainWindow>

#include "cartocrow/core/core.h"
#include "cartocrow/core/timer.h"
#include "cartocrow/stenomap/medialAxis.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::medial_axis;

StenomapDemo::StenomapDemo() {
	setWindowTitle("AAACartoCrow – Stenomap demo");

  // Make simple polygon
  Polygon<Inexact> polygon;
	polygon.push_back(Point<Inexact>(80, 50));
	polygon.push_back(Point<Inexact>(90, 4));
	polygon.push_back(Point<Inexact>(40, 30));
	polygon.push_back(Point<Inexact>(60, 60));
  m_polygons.push_back(polygon);

  // setup renderer
	m_renderer = new GeometryWidget();
	m_renderer->setMaxZoom(10);
	m_renderer->setGridMode(GeometryWidget::GridMode::CARTESIAN);
	setCentralWidget(m_renderer);

	m_medialAxisBox = new QCheckBox("Compute with obstacles");
	connect(m_medialAxisBox, &QCheckBox::stateChanged, [&]() {
		recalculate();
	});
	QToolBar* toolBar = new QToolBar();
	toolBar->addWidget(m_medialAxisBox);

  recalculate();
}

void StenomapDemo::recalculate() {
  // draw polygon
  // TODO: figure out how to draw a simple polygon
  GeometryRenderer rend;
  for (const Polygon<Inexact>& p : m_polygons) {
    rend->draw(p);
    /* GeometryRenderer::draw(p); */
  }

	if (m_medialAxisBox->isChecked()) {
    // find/compute medial axis and draw it
    for (const Polygon<Inexact>& p : m_polygons) {
      MedialAxis medial_axis = MedialAxis(p);
      rend->draw(medial_axis);
      /* GeometryRenderer::draw(p); */
    }
  }

}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	StenomapDemo demo;
	demo.show();
	app.exec();
}
