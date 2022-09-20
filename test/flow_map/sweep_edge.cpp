#include "../catch.hpp"

#include "../../cartocrow/flow_map/sweep_edge.h"

using namespace cartocrow;
using namespace cartocrow::flow_map;

TEST_CASE("Intersecting a segment and a left spiral") {
	const SweepEdgeShape segment(PolarPoint(2, M_PI / 2), PolarPoint(Point<Inexact>(-6, 8)));
	const SweepEdgeShape spiral(SweepEdgeShape::Type::LEFT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	auto r = segment.intersectWith(spiral, 2);
	REQUIRE(r.has_value());
	CHECK(spiral.phiForR(*r) == Approx(segment.phiForR(*r)));
}

TEST_CASE("Intersecting a segment and a left spiral without intersections") {
	const SweepEdgeShape segment(PolarPoint(2, M_PI / 2), PolarPoint(Point<Inexact>(6, 8)));
	const SweepEdgeShape spiral(SweepEdgeShape::Type::LEFT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	auto r = segment.intersectWith(spiral, 2);
	CHECK(!r.has_value());
}

TEST_CASE("Intersecting a segment and a right spiral") {
	const SweepEdgeShape segment(PolarPoint(2, M_PI / 2), PolarPoint(Point<Inexact>(6, 8)));
	const SweepEdgeShape spiral(SweepEdgeShape::Type::RIGHT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	auto r = segment.intersectWith(spiral, 2);
	REQUIRE(r.has_value());
	CHECK(spiral.phiForR(*r) == Approx(segment.phiForR(*r)));
}

TEST_CASE("Intersecting a segment and a right spiral without intersections") {
	const SweepEdgeShape segment(PolarPoint(2, M_PI / 2), PolarPoint(Point<Inexact>(-6, 8)));
	const SweepEdgeShape spiral(SweepEdgeShape::Type::RIGHT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	auto r = segment.intersectWith(spiral, 2);
	CHECK(!r.has_value());
}

TEST_CASE("Intersecting a left spiral and a right spiral starting at the same point") {
	const SweepEdgeShape spiral1(SweepEdgeShape::Type::RIGHT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	const SweepEdgeShape spiral2(SweepEdgeShape::Type::LEFT_SPIRAL, PolarPoint(2, M_PI / 2), 0.5);
	auto r = spiral1.intersectWith(spiral2, 2);
	REQUIRE(r.has_value());
	CHECK(r > 2);
	CHECK(spiral1.phiForR(*r) == Approx(spiral2.phiForR(*r)));
}

TEST_CASE("Intersecting a left spiral and a right spiral starting at different points") {
	const SweepEdgeShape spiral1(SweepEdgeShape::Type::RIGHT_SPIRAL, PolarPoint(2, 3 * M_PI / 4),
	                             0.5);
	const SweepEdgeShape spiral2(SweepEdgeShape::Type::LEFT_SPIRAL, PolarPoint(2, M_PI / 4), 0.5);
	auto r = spiral1.intersectWith(spiral2, 2);
	REQUIRE(r.has_value());
	CHECK(r > 2);
	CHECK(spiral1.phiForR(*r) == Approx(spiral2.phiForR(*r)));
}