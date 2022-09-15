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

#include "spiral_tree_demo.h"

#include <QApplication>

#include "cartocrow/core/core.h"
#include "cartocrow/flow_map/painting.h"
#include "cartocrow/flow_map/parameters.h"
#include "cartocrow/flow_map/place.h"
#include "cartocrow/flow_map/spiral_tree.h"
#include "cartocrow/flow_map/spiral_tree_obstructed_algorithm.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/renderer/ipe_renderer.h"

using namespace cartocrow;
using namespace cartocrow::flow_map;
using namespace cartocrow::renderer;

SpiralTreeDemo::SpiralTreeDemo() {
	setWindowTitle("CartoCrow – Spiral tree demo");

	/*
	m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(8, 95));
	m_obstacle.push_back(Point<Inexact>(50, 140));
	m_obstacle.push_back(Point<Inexact>(-43, 134));
	m_obstacle.push_back(Point<Inexact>(-50, 100));
	*/

	m_obstacle.push_back(Point<Inexact>(-3, 38));
	m_obstacle.push_back(Point<Inexact>(11, 56));
	m_obstacle.push_back(Point<Inexact>(21, 79));
	m_obstacle.push_back(Point<Inexact>(15, 108));
	m_obstacle.push_back(Point<Inexact>(-57, 140));

	// join not working
	/*m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(20, 80));
	m_obstacle.push_back(Point<Inexact>(-20, 70));*/

	m_renderer = new GeometryWidget();
	setCentralWidget(m_renderer);

	recalculate();
	connect(m_renderer, &GeometryWidget::dragStarted, [&](Point<Inexact> p) {
		m_draggedPoint = findClosestPoint(p, 10 / m_renderer->zoomFactor());
		recalculate();
	});
	connect(m_renderer, &GeometryWidget::dragMoved, [&](Point<Inexact> p) {
		if (m_draggedPoint != nullptr) {
			Point<Inexact> originalPoint = *m_draggedPoint;
			*m_draggedPoint = p;
			if (!m_obstacle.is_simple()) {
				*m_draggedPoint = originalPoint;
			}
			recalculate();
		}
	});
	connect(m_renderer, &GeometryWidget::dragEnded, [&](Point<Inexact> p) {
		m_draggedPoint = nullptr;
		recalculate();
	});
}

void SpiralTreeDemo::recalculate() {
	auto tree = std::make_shared<SpiralTree>(Point<Inexact>(0, 0), 0.5061454830783556);
	tree->addPlace("p1", Point<Inexact>(0, 400), 1);
	tree->addObstacle(m_obstacle);
	SpiralTreeObstructedAlgorithm algorithm(tree);
	algorithm.run();

	Painting::Options options;
	Painting p(nullptr, tree, options);
	auto painting = std::make_shared<Painting>(nullptr, tree, options);

	m_renderer->clear();
	m_renderer->addPainting(painting);
	m_renderer->addPainting(algorithm.debugPainting());

	m_renderer->update();
}

Point<Inexact>* SpiralTreeDemo::findClosestPoint(Point<Inexact> p, Number<Inexact> radius) {
	Point<Inexact>* closest = nullptr;
	Number<Inexact> minSquaredDistance = radius * radius;
	for (auto& vertex : m_obstacle) {
		Number<Inexact> squaredDistance = (vertex - p).squared_length();
		if (squaredDistance < minSquaredDistance) {
			minSquaredDistance = squaredDistance;
			closest = &vertex;
		}
	}
	return closest;
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	SpiralTreeDemo demo;
	demo.show();
	app.exec();
}