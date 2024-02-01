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

	/*m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(25, 40));
	m_obstacle.push_back(Point<Inexact>(8, 95));
	m_obstacle.push_back(Point<Inexact>(50, 140));
	m_obstacle.push_back(Point<Inexact>(-41, 134));
	m_obstacle.push_back(Point<Inexact>(-43, 60));
	m_obstacle.push_back(Point<Inexact>(0, 70));
	m_obstacle.push_back(Point<Inexact>(-20, 20));*/

	// join not working
	/*m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(20, 80));
	m_obstacle.push_back(Point<Inexact>(-20, 70));*/

	// right vertex event bug
	/*m_obstacle.push_back(Point<Inexact>(-3, 38));
	m_obstacle.push_back(Point<Inexact>(11, 56));
	m_obstacle.push_back(Point<Inexact>(21, 79));
	m_obstacle.push_back(Point<Inexact>(15, 108));
	m_obstacle.push_back(Point<Inexact>(-57, 140));*/

	// far vertex event on near side bug
	/*m_obstacle.push_back(Point<Inexact>(-19, 59));
	m_obstacle.push_back(Point<Inexact>(9, 83));
	m_obstacle.push_back(Point<Inexact>(26, 59));
	m_obstacle.push_back(Point<Inexact>(49, 114));
	m_obstacle.push_back(Point<Inexact>(-33, 117));*/

	// another far vertex event handling bug
	/*m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(8, 95));
	m_obstacle.push_back(Point<Inexact>(50, 140));
	m_obstacle.push_back(Point<Inexact>(0, 175));
	m_obstacle.push_back(Point<Inexact>(-50, 100));*/

	// join event with obstacle edge
	/*m_obstacle.push_back(Point<Inexact>(0, 50));
	m_obstacle.push_back(Point<Inexact>(8, 95));
	m_obstacle.push_back(Point<Inexact>(50, 140));
	m_obstacle.push_back(Point<Inexact>(-37, 175));
	m_obstacle.push_back(Point<Inexact>(-50, 100));*/

	// φ = π issues
	/*m_obstacle.push_back(Point<Inexact>(-7, 19));
	m_obstacle.push_back(Point<Inexact>(-76, 30));
	m_obstacle.push_back(Point<Inexact>(-82, -47));*/

	// hitting the same obstacle edge with two spirals
	/*m_obstacle.push_back(Point<Inexact>(0, 7));
	m_obstacle.push_back(Point<Inexact>(4.5, 11));
	m_obstacle.push_back(Point<Inexact>(-5.5, 9));*/

	m_places.push_back(std::make_shared<Point<Inexact>>(11.2121212, 17.0707070));
	m_places.push_back(std::make_shared<Point<Inexact>>(13.9393939, -14.1414141));
	m_places.push_back(std::make_shared<Point<Inexact>>(-4.5454545, -18.9898989));
	m_places.push_back(std::make_shared<Point<Inexact>>(16.6666666, 6.1616161));
	m_places.push_back(std::make_shared<Point<Inexact>>(-9.8989898, 13.9393939));
	m_places.push_back(std::make_shared<Point<Inexact>>(-16.1616161, -2.6262626));

	auto obstacle = std::make_shared<Polygon<Inexact>>();
	obstacle->push_back(Point<Inexact>(-8.9898989, -5.4545454));
	obstacle->push_back(Point<Inexact>(-9.5959595, -0.4040404));
	obstacle->push_back(Point<Inexact>(-4.7474747, -3.0303030));
	obstacle->push_back(Point<Inexact>(-6.5656565, -6.7676767));
	m_obstacles.push_back(obstacle);

	obstacle = std::make_shared<Polygon<Inexact>>();
	obstacle->push_back(Point<Inexact>(1.3131313, 10.2020202));
	obstacle->push_back(Point<Inexact>(6.1616161, 10.4040404));
	obstacle->push_back(Point<Inexact>(5.6565656, 5.2525252));
	m_obstacles.push_back(obstacle);

	obstacle = std::make_shared<Polygon<Inexact>>();
	obstacle->push_back(Point<Inexact>(4.6464646, -10.4040404));
	obstacle->push_back(Point<Inexact>(10.4040404, -7.1717171));
	obstacle->push_back(Point<Inexact>(7.4747474, -13.9393939));
	m_obstacles.push_back(obstacle);

	m_renderer = new GeometryWidget();
	m_renderer->setMaxZoom(10000);
	m_renderer->setGridMode(GeometryWidget::GridMode::POLAR);
	setCentralWidget(m_renderer);

	QToolBar* toolBar = new QToolBar();
	m_obstacleBox = new QCheckBox("Compute with obstacles");
	connect(m_obstacleBox, &QCheckBox::stateChanged, [&]() {
		recalculate();
	});
	toolBar->addSeparator();
	toolBar->addWidget(m_obstacleBox);
	toolBar->addWidget(new QLabel("α = "));
	m_alphaSlider = new QSlider(Qt::Horizontal);
	m_alphaSlider->setMinimum(0);
	m_alphaSlider->setMaximum(499);
	m_alphaSlider->setValue(139);
	toolBar->addWidget(m_alphaSlider);
	addToolBar(toolBar);
	m_alphaLabel = new QLabel("0.139π");
	toolBar->addWidget(m_alphaLabel);
	connect(m_alphaSlider, &QSlider::valueChanged, [&](int value) {
		m_alpha = M_PI * value / 1000.0;
		m_alphaLabel->setText(QString("%1π").arg(value / 1000.0, 0, 'f', 3));
		recalculate();
	});
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
	for (const auto& place : m_places) {
		tree->addPlace("", *place, 1);
	}
	for (auto obstacle : m_obstacles) {
		tree->addObstacle(*obstacle);
	}
	tree->addShields();
	t.stamp("Constructing tree and obstacles");

	m_renderer->clear();
	if (m_obstacleBox->isChecked()) {
		ReachableRegionAlgorithm reachableRegionAlg(tree);
		ReachableRegionAlgorithm::ReachableRegion reachableRegion = reachableRegionAlg.run();
		t.stamp("Computing reachable region");

		SpiralTreeObstructedAlgorithm spiralTreeAlg(tree, reachableRegion);
		spiralTreeAlg.run();
		t.stamp("Computing spiral tree");

		m_renderer->addPainting(reachableRegionAlg.debugPainting(), "Reachable region sweep");
		m_renderer->addPainting(spiralTreeAlg.debugPainting(), "Spiral tree sweep");
	} else {
		SpiralTreeUnobstructedAlgorithm spiralTreeAlg(*tree);
		spiralTreeAlg.run();
		t.stamp("Computing spiral tree");

		m_renderer->addPainting(spiralTreeAlg.debugPainting(), "Spiral tree sweep");
	}

	t.output();

	Painting::Options options;
	auto painting = std::make_shared<Painting>(nullptr, tree, options);
	m_renderer->addPainting(painting, "Spiral tree");

	m_renderer->update();
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	StenomapDemo demo;
	demo.show();
	app.exec();
}
