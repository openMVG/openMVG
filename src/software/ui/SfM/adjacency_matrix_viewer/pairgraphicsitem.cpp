// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "pairgraphicsitem.h"

#include <QtWidgets>

PairGraphicsItem::PairGraphicsItem
(
  const QColor &color,
  unsigned int x,
  unsigned int y,
  unsigned int matches_count
)
: color(color),
  x(x),
  y(y),
  matches_count(matches_count)
{
    setFlags(ItemIsSelectable);
    setAcceptHoverEvents(true);
}

QRectF PairGraphicsItem::boundingRect() const
{
    return QRectF(0, 0, 100, 100);
}

QPainterPath PairGraphicsItem::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void PairGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
  Q_UNUSED(widget);

  QColor fillColor = (option->state & QStyle::State_Selected) ? color.dark(150) : color;
  if (option->state & QStyle::State_MouseOver)
    fillColor = fillColor.light(125);

  const qreal lod = option->levelOfDetailFromTransform(painter->worldTransform());

  if (lod < 0.5)
  {
    painter->fillRect(QRectF(0, 0, 100, 100), fillColor);
    return;
  }
  else
  {
    const QBrush b = painter->brush();

    painter->setBrush(QBrush(fillColor.dark(option->state & QStyle::State_Sunken ? 120 : 100)));
    painter->drawRect(QRect(0, 0, 100, 100));
    painter->setBrush(b);

    // Draw text (information about the camera pair)
    if (lod < 2)
    {
      QFont font("Times", 10);
      font.setStyleStrategy(QFont::ForceOutline);
      painter->setFont(font);
      painter->save();
      painter->scale(1.2, 1.2);
      painter->drawText(10, 20, QString("%1").arg(x));
      painter->drawText(10, 60, QString("%1").arg(y));
      painter->restore();
    }
    else // (lod >= 2)
    {
        QFont font("Times", 10);
        font.setStyleStrategy(QFont::ForceOutline);
        painter->setFont(font);
        painter->save();
        painter->scale(0.8, 0.8);
        painter->drawText(10, 20, QString("Pair:     %1; %2").arg(x).arg(y));
        painter->drawText(10, 60, QString("#Matches: %1").arg(matches_count));
        painter->restore();
    }
  }
}
