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
#include "grid_painting.h"

#include <QApplication>
#include <QCheckBox>
#include <QMainWindow>

#include "cartocrow/core/core.h"
#include "cartocrow/core/timer.h"
#include "cartocrow/stenomap/medial_axis.h"
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
    // ma.calculate_weight_function();
    /* MedialAxis ma(polygon);
       ma.print_adjacency_list();
       ma.compute_branches();
       ma.remove_branch(0);
       ma.print_adjacency_list(); */
}

void StenomapDemo::recalculate() {
    for (const Polygon<Inexact>& p : m_polygons) {
        // Draw polygon
        m_renderer->addPainting(std::make_shared<PolygonPainting>(PolygonPainting(p)), "Polygon");

        MedialAxis medial_axis(p); 
        medial_axis.compute_grid(10, 10, 5);

        /* auto medialAxisPainting = std::make_shared<MedialAxisPainting>(medial_axis);
        auto grid_painting = std::make_shared<GridPainting>(medial_axis); */
        m_renderer->addPainting(std::make_shared<MedialAxisPainting>(medial_axis), "Medial Axis");
        m_renderer->addPainting(std::make_shared<GridPainting>(medial_axis), "Grid");
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    StenomapDemo demo;
    demo.show();
    app.exec();
}
