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

#include <QApplication>
#include <QCheckBox>
#include <QMainWindow>

#include "cartocrow/core/core.h"
#include "cartocrow/core/timer.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/renderer/geometry_renderer.h"

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::medial_axis;

StenomapDemo::StenomapDemo() {
    setWindowTitle("CartoCrow – Stenomap demo");

    // Make simple polygon
    Polygon<Inexact> polygon;
    /* polygon.push_back(Point<Inexact>(80, 50));
       polygon.push_back(Point<Inexact>(90, 4));
       polygon.push_back(Point<Inexact>(40, 30));
       polygon.push_back(Point<Inexact>(60, 60));*/ 
    polygon.push_back( Point<Inexact>(-1,-1) ) ;
    polygon.push_back( Point<Inexact>(0,-12) ) ;
    polygon.push_back( Point<Inexact>(1,-1) ) ;
    polygon.push_back( Point<Inexact>(12,0) ) ;
    polygon.push_back( Point<Inexact>(1,1) ) ;
    polygon.push_back( Point<Inexact>(0,12) ) ;
    polygon.push_back( Point<Inexact>(-1,1) ) ;
    polygon.push_back( Point<Inexact>(-12,0) ) ;
    m_polygons.push_back(polygon);

    MedialAxis medial_axis(polygon);
    medial_axis.compute_grid(10, 10, 5);
    // medial_axis.prune_grid();
    // medial_axis.calculate_weight_function();
    // medial_axis.compute_centroid_neighborhoods();
    // medial_axis.compute_centroid_closest_points();
    medial_axis.compute_branches();
    medial_axis.compute_grid_closest_branches();
    medial_axis.print_grid_closest_branches();

    MedialAxis old = medial_axis;
    /* medial_axis.prune_points(0.05); */

    // setup renderer
    m_renderer = new GeometryWidget();
    m_renderer->setMaxZoom(10);
    m_renderer->setGridMode(GeometryWidget::GridMode::CARTESIAN);
    setCentralWidget(m_renderer);

    m_medialAxisBox = new QCheckBox("Compute with obstacles");
    connect(m_medialAxisBox, &QCheckBox::stateChanged, [&]() {
            recalculate(medial_axis, old);
            });
    QToolBar* toolBar = new QToolBar();
    toolBar->addWidget(m_medialAxisBox);

    recalculate(medial_axis, old);
}

void StenomapDemo::recalculate(const MedialAxis medial_axis, const MedialAxis old) {
    for (const Polygon<Inexact>& p : m_polygons) {
        // Draw polygon
        m_renderer->addPainting(std::make_shared<PolygonPainting>(PolygonPainting(p)), "Polygon");
        // Draw medial axis
        m_renderer->addPainting(std::make_shared<MedialAxisPainting>(old), "Medial Axis");
        // Draw pruned medial axis
        m_renderer->addPainting(std::make_shared<PrunedMedialAxisPainting>(medial_axis), "Pruned medial Axis");
        // Draw grid
        m_renderer->addPainting(std::make_shared<GridPainting>(medial_axis), "Grid");
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    StenomapDemo demo;
    demo.show();
    app.exec();
}
