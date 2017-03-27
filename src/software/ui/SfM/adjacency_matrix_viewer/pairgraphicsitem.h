// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <QColor>
#include <QGraphicsItem>

//
// QTGraphicsItem to represent an image pair
//
class PairGraphicsItem : public QGraphicsItem
{
public:
    PairGraphicsItem
    (
      const QColor &color,
      unsigned int x, // View index
      unsigned int y, // View index
      unsigned int matches_count
    );

    void paint
    (
      QPainter *painter,
      const QStyleOptionGraphicsItem *item,
      QWidget *widget
    ) Q_DECL_OVERRIDE;

    QRectF boundingRect() const Q_DECL_OVERRIDE;
    QPainterPath shape() const Q_DECL_OVERRIDE;

    unsigned int get_x() const {return x;}
    unsigned int get_y() const {return y;}

private:
    unsigned int x;
    unsigned int y;
    unsigned int matches_count;
    QColor color;
};
