// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ControlPoint2DNode.hpp"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <iostream>

ControlPoint2DNode::ControlPoint2DNode
(
  const QPointF& pos,
  double & x,
  double & y,
  size_t id_control_point
):
  _id_control_point(id_control_point),
  _x(x),
  _y(y)
{
  setFlags(ItemIsMovable | ItemIsSelectable);
  setFlag(ItemSendsGeometryChanges);
  setCacheMode(DeviceCoordinateCache);
  setZValue(1.);
  setPos(pos);
}

size_t ControlPoint2DNode::controlPointId() const
{
  return _id_control_point;
}

QRectF ControlPoint2DNode::boundingRect() const
{
  return QRectF(QPointF(-10,-10),QSize(20,20));
}

QPainterPath ControlPoint2DNode::shape() const
{
  QPainterPath path;
  path.addEllipse(QPointF(0,0),10,10);
  return path;
}

void ControlPoint2DNode::paint
(
  QPainter *painter,
  const QStyleOptionGraphicsItem *option,
  QWidget *widget
)
{
  QRadialGradient gradient(-3, -3, 10);
  if (isSelected())  {
    gradient.setCenter(3, 3);
    gradient.setFocalPoint(3, 3);
    gradient.setColorAt(1, QColor(Qt::yellow).light(120));
    gradient.setColorAt(0, QColor(Qt::darkYellow).light(120));
  }
  else  {
    gradient.setColorAt(0, Qt::yellow);
    gradient.setColorAt(1, Qt::darkYellow);
  }
  painter->setBrush(gradient);
  painter->setPen(QPen(Qt::black, 0));
  painter->drawEllipse(QPointF(0,0),10,10);

  // Or just draw a cross ?
  //painter->setPen(QPen(Qt::black, 1));
  //painter->drawLine(10, 10, -10, -10);
  //painter->drawLine(10, -10, -10, 10);

  // Paint the index of the control_point
  // adapt the font size to fit the GCP id to the bounding box of the ellipse
  const QRectF node_rect(QPointF(-8,-8),QSize(18,18));
  const QString s_node_id = QString::number(_id_control_point);
  const float factor = node_rect.width() / painter->fontMetrics().width(s_node_id);
  painter->setPen(QPen(Qt::black, 5));
  if (factor < 1)
  {
    QFont f = painter->font();
    f.setPointSizeF(f.pointSizeF()*factor);
    painter->setFont(f);
  }
  painter->drawText(node_rect, Qt::AlignCenter, s_node_id);
}

QVariant ControlPoint2DNode::itemChange
(
  GraphicsItemChange change,
  const QVariant &value
)
{
  const QVariant variant = QGraphicsItem::itemChange(change, value);
  _x = scenePos().x();
  _y = scenePos().y();
  return variant;
}

void ControlPoint2DNode::mousePressEvent
(
  QGraphicsSceneMouseEvent *event
)
{
  update();
  QGraphicsItem::mousePressEvent(event);
}
void ControlPoint2DNode::mouseReleaseEvent
(
  QGraphicsSceneMouseEvent *event
)
{
  update();
  QGraphicsItem::mouseReleaseEvent(event);
}
