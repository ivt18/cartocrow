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

#include "medial_axis_demo.h"
#include "polygon_painting.h"
#include "medial_axis_painting.h"
#include "pruned_medial_axis_painting.h"
#include "grid_painting.h"
#include "pruned_grid_painting.h"

#include <QApplication>
#include <QCheckBox>
#include <QMainWindow>

#include "cartocrow/core/core.h"
#include "cartocrow/core/timer.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/renderer/geometry_renderer.h"

#include "cartocrow/core/ipe_reader.h"
#include "cartocrow/core/region_map.h"
#include <ipedoc.h>
#include <ipepath.h>

#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::medial_axis;

StenomapDemo::StenomapDemo() {
    setWindowTitle("CartoCrow – Stenomap demo");
    std::filesystem::path ipePath = "/home/diego/src/cartocrow/data/europe.ipe";
	RegionMap map = ipeToRegionMap(ipePath);
    for (const auto a : map) {
        std::cout << a.first << std::endl;
    }
    CGAL::Polyline_simplification_2::Squared_distance_cost cost;
    CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold stop(0.1);
    Polygon<Inexact> pgon;
    PolygonSet<Exact> polygonSet;
	std::vector<PolygonWithHoles<Exact>> polygons;

    polygonSet = map["ESP"].shape;
    polygons.clear();
	polygonSet.polygons_with_holes(std::back_inserter(polygons));
    double maxArea = 0;
    // To get mainland Spain we get the biggest polygon from the set (instead of some random island)
	for (const PolygonWithHoles<Exact> p : polygons) {
        if (approximate(p.outer_boundary()).area() > maxArea) {
            pgon = approximate(p.outer_boundary());
            maxArea = pgon.area();
        }
	}

    // Simplifying polygon
    pgon = CGAL::Polyline_simplification_2::simplify(pgon, cost, stop);
    m_polygons.push_back(pgon);

    polygonSet = map["FRA"].shape;
	polygons.clear();
	polygonSet.polygons_with_holes(std::back_inserter(polygons));
    maxArea = 0;
    // To get mainland France we get the biggest polygon from the set
	for (const PolygonWithHoles<Exact> p : polygons) {
        if (approximate(p.outer_boundary()).area() > maxArea) {
            pgon = approximate(p.outer_boundary());
            maxArea = pgon.area();
        }
	}
    // Simplifying polygon
    CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold stop2(0.5);
    pgon = CGAL::Polyline_simplification_2::simplify(pgon, cost, stop);
    m_polygons.push_back(pgon);

    /* Polygon<Inexact> polygon; */
    /* // Make simple polygon */
    /* polygon.push_back( Point<Inexact>(-1,-1) ) ; */
    /* polygon.push_back( Point<Inexact>(-4,-22) ) ; */
    /* polygon.push_back( Point<Inexact>(0,-22) ) ; */
    /* polygon.push_back( Point<Inexact>(1,-1) ) ; */
    /* polygon.push_back( Point<Inexact>(32,0) ) ; */
    /* polygon.push_back( Point<Inexact>(32,20) ) ; */
    /* polygon.push_back( Point<Inexact>(1,1) ) ; */
    /* polygon.push_back( Point<Inexact>(0,42) ) ; */
    /* polygon.push_back( Point<Inexact>(-1,1) ) ; */
    /* polygon.push_back( Point<Inexact>(-12,0) ) ; */
    /* m_polygons.push_back(polygon); */

    for (Polygon<Inexact> polygon : m_polygons) {
        MedialAxis medial_axis(polygon);

        medial_axis.compute_grid(50, 50);
        medial_axis.prune_grid();
        medial_axis.compute_branches();

        m_medialAxisOld.push_back(medial_axis);
        medial_axis.prune_points(0.015);
        m_medialAxis.push_back(medial_axis);
    }
    // setup renderer
    m_renderer = new GeometryWidget();
    m_renderer->setMaxZoom(1000);
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
    for (int i = 0; i < m_polygons.size(); i++) {
        std::string index = std::to_string(i);
        // Draw polygon
        m_renderer->addPainting(std::make_shared<PolygonPainting>(PolygonPainting(m_polygons[i])), "Polygon " + index);
        // Draw medial axis
        m_renderer->addPainting(std::make_shared<MedialAxisPainting>(m_medialAxisOld[i]), "Medial Axis" + index);
        // Draw pruned medial axis
        m_renderer->addPainting(std::make_shared<PrunedMedialAxisPainting>(m_medialAxis[i]), "Pruned medial Axis" + index);
        // Draw grid
        m_renderer->addPainting(std::make_shared<GridPainting>(m_medialAxis[i]), "Grid" + index);
        // Draw pruned grid
        m_renderer->addPainting(std::make_shared<PrunedGridPainting>(m_medialAxis[i]), "Pruned grid" + index);
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    StenomapDemo demo;
    demo.show();
    app.exec();
}
