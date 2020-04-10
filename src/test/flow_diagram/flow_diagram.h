/*
The Necklace Map console application implements the algorithmic
geo-visualization method by the same name, developed by
Bettina Speckmann and Kevin Verbeek at TU Eindhoven
(DOI: 10.1109/TVCG.2010.180 & 10.1142/S021819591550003X).
Copyright (C) 2019  Netherlands eScience Center and TU Eindhoven

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

Created by tvl (t.vanlankveld@esciencecenter.nl) on 10-04-2020
*/

#ifndef GEOVIZ_TEST_FLOW_DIAGRAM_FLOW_DIAGRAM_H
#define GEOVIZ_TEST_FLOW_DIAGRAM_FLOW_DIAGRAM_H

#include "test/test.h"


void TestFlowDiagram() {}  // Linking hack, each new test cpp file has it.

SUITE(suite_FlowDiagram) {

TEST(SanityFlowDiagram)
{
  CHECK_EQUAL(2, 2);
}

} // suite_FlowDiagram

#endif //GEOVIZ_TEST_FLOW_DIAGRAM_FLOW_DIAGRAM_H