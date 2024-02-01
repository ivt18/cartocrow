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
#include "cartocrow/flow_map/painting.h"
#include "cartocrow/flow_map/parameters.h"
#include "cartocrow/flow_map/place.h"
#include "cartocrow/flow_map/reachable_region_algorithm.h"
#include "cartocrow/flow_map/spiral_tree.h"
#include "cartocrow/flow_map/spiral_tree_obstructed_algorithm.h"
#include "cartocrow/flow_map/spiral_tree_unobstructed_algorithm.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"

using namespace cartocrow;
using namespace cartocrow::flow_map;
using namespace cartocrow::renderer;

StenomapDemo::StenomapDemo() {
	setWindowTitle("CartoCrow – Stenomap demo");

	auto obstacle = std::make_shared<Polygon<Inexact>>();
	obstacle->push_back(Point<Inexact>(-80.9898989, 50.4545454));
	obstacle->push_back(Point<Inexact>(-90.5959595, 4.040404));
	obstacle->push_back(Point<Inexact>(-40.7474747, 30.0303030));
	obstacle->push_back(Point<Inexact>(-60.5656565, 60.7676767));
	m_obstacles.push_back(obstacle);

	m_renderer = new GeometryWidget();
	m_renderer->setMaxZoom(10);
	m_renderer->setGridMode(GeometryWidget::GridMode::CARTESIAN);
	setCentralWidget(m_renderer);

	for (auto place : m_places) {
		m_renderer->registerEditable(place);
	}
	for (auto obstacle : m_obstacles) {
		m_renderer->registerEditable(obstacle);
	}
	connect(m_renderer, &GeometryWidget::edited, [&]() {
		recalculate();
	});
	recalculate();
}

void StenomapDemo::recalculate() {
	Timer t;
  auto tree = std::make_shared<SpiralTree>(Point<Inexact>(0, 0), m_alpha);
	t.stamp("Constructing polygon");

  for (auto obstacle : m_obstacles) {
		tree->addObstacle(*obstacle);
	}

	m_renderer->clear();

	t.output();

	Painting::Options options;
	auto painting = std::make_shared<Painting>(nullptr, tree, options);
	m_renderer->addPainting(painting, "Polygon");

	m_renderer->update();
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	StenomapDemo demo;
	demo.show();
	app.exec();
}
