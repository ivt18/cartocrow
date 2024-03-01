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

  MedialAxis ma(polygon);
  ma.compute_grid(10, 10, 5);
  ma.calculate_weight_function();
  ma.print_adjacency_list();
    /* MedialAxis ma(polygon);
    ma.print_adjacency_list();
    ma.compute_branches();
    ma.remove_branch(0);
    ma.print_adjacency_list(); */
}

void StenomapDemo::recalculate() {
    // draw polygon
    for (const Polygon<Inexact>& p : m_polygons) {
        m_renderer->addPainting(std::make_shared<PolygonPainting>(PolygonPainting(p)), "Polygon");
        std::cout<<"ok";
    }

    // TODO: make sure this works well with medial axis computation implementation once 
    // drawing skeleton is added and feature/feature-points is merged

   // if (m_medialAxisBox->isChecked()) {
        // find/compute medial axis and draw it
        for (const Polygon<Inexact>& p : m_polygons) {
            std::cout<<"ok1";
            // TODO: get SsPtr medialAxis from p using medial_axis.h
            /* MedialAxis med_axis = MedialAxis(p); */
            /* MedialAxisPainting m_painting = MedialAxisPainting(med_axis.iss); */
            /* m_renderer->addPainting(std::make_shared<PolygonPainting>(m_painting), "medialAxis"); */
         MedialAxis med_axis(p); 
        auto medialAxisPainting = std::make_shared<MedialAxisPainting>(med_axis.getIss());
        m_renderer->addPainting(medialAxisPainting, "Medial Axis");
        }
   // }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    StenomapDemo demo;
    demo.show();
    app.exec();
}
