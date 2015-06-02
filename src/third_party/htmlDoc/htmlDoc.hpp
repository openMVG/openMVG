
/**
 * @file htmlDoc.h
 * @brief Simple html document writer
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef MIMATTE_HTML_DOC_H
#define MIMATTE_HTML_DOC_H

#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <vector>

namespace htmlDocument
{
  inline const std::string htmlMarkup(const std::string & markup, const std::string & text)
  {
    std::ostringstream os;
    os << '<'<< markup <<'>' << text << "</"<< markup <<'>' <<"\n";
    return os.str();
  }

  inline const std::string htmlComment(const std::string & text)
  {
    std::ostringstream os;
    os << "<!-- "<<  text << " -->" << "\n";
    return os.str();
  }

  inline const std::string htmlOpenMarkup(const std::string & markup,
    const std::string & attributes)
  {
    std::ostringstream os;
    os << '<' << markup << ' ' << attributes  << "/>"  <<"\n";
    return os.str();
  }

  /// Return a chain in the form attributes="val"
  template<typename T>
  inline const std::string quotedAttributes(const std::string & attributes, const T & val)
  {
    std::ostringstream os;
    os << attributes << "=\"" << val << '"';
    return os.str();
  }

  /// Return a chain of the value T
  template<typename T>
  inline const std::string toString(const T & val)
  {
    std::ostringstream os;
    os << val;
    return os.str();
  }

  template<typename T, typename T2>
  static std::pair< std::pair<T,T>, std::pair<T,T> > autoJSXGraphViewport(const std::vector<T> & vec_x, const std::vector<T2> & vec_y)
  {
    std::pair< std::pair<T,T>, std::pair<T,T> > viewport = std::make_pair( std::make_pair(0,0), std::make_pair(0,0));
    if (!vec_x.empty() && !vec_y.empty() && vec_x.size() == vec_y.size())
    {
      T minValX, maxValX;
      //For X values
      minValX = *min_element(vec_x.begin(), vec_x.end());
      maxValX = *max_element(vec_x.begin(), vec_x.end());

      //For Y values
      T2 minValY, maxValY;
      minValY = *min_element(vec_y.begin(), vec_y.end());
      maxValY = *max_element(vec_y.begin(), vec_y.end());

      //Use the value with a little margin
      T rangeX = maxValX-minValX;
      T rangeY = maxValY-minValY;
      viewport = std::make_pair( std::make_pair
        (-.2*rangeX+minValX,0.2*rangeX+maxValX),
        std::make_pair(-.2*rangeY+minValY,0.2*rangeY+maxValY));
    }
    return viewport;
  }

  template<typename T, typename T2>
  static std::pair< std::pair<T,T>, std::pair<T,T> > autoJSXGraphViewport(const std::vector<T2> & vec_y, bool bForceY0 = true)
  {
    std::pair< std::pair<T,T>, std::pair<T,T> > viewport = std::make_pair( std::make_pair(0,0), std::make_pair(0,0));
    if (!vec_y.empty())
    {
      T2 minValX, maxValX;
      //For X values
      minValX = 0;
      maxValX = static_cast<T2>(vec_y.size());

      //For Y values
      T2 minValY, maxValY;
      minValY = *min_element(vec_y.begin(), vec_y.end());
      maxValY = *max_element(vec_y.begin(), vec_y.end());

      if (bForceY0)  {  minValY = T2(0);  }
      //Use the value with a little margin
      T2 rangeX = maxValX-minValX;
      T2 rangeY = maxValY-minValY;
      viewport = std::make_pair( std::make_pair
        (static_cast<T>(-.2*rangeX+minValX),static_cast<T>(0.2*rangeX+maxValX)),
        std::make_pair(static_cast<T>(-.2*rangeY+minValY),static_cast<T>(0.2*rangeY+maxValY)));
    }
    return viewport;
  }

  class htmlDocumentStream
  {
  public:
    htmlDocumentStream(const std::string & title)
    {
      htmlStream << htmlMarkup("head",
        std::string("\n") +
        "<link rel=\"stylesheet\" type=\"text/css\" href=\"http://jsxgraph.uni-bayreuth.de/distrib/jsxgraph.css\" />\n" +
        "<link rel=\"stylesheet\" type=\"text/css\" href=\"http://imagine.enpc.fr/~moulonp/style.css\" />" +
        "<script type=\"text/javascript\" src=\"http://cdnjs.cloudflare.com/ajax/libs/jsxgraph/0.93/jsxgraphcore.js\"></script>\n" +
        htmlMarkup("title",title));
    }

    htmlDocumentStream(const std::string & title,
                       const std::vector<std::string> & vec_css,
                       const std::vector<std::string> & vec_js)
    {
      htmlStream << "\n<head>\n";
      htmlStream << htmlMarkup("title",title);
      // CSS and JS ressources
      for (std::vector<std::string>::const_iterator iter = vec_css.begin(); iter != vec_css.end(); ++iter)
        htmlStream << "<link rel=\"stylesheet\" type=\"text/css\" href=\"" << *iter <<"\" />\n";
      for (std::vector<std::string>::const_iterator iter = vec_js.begin(); iter != vec_js.end(); ++iter)
        htmlStream << "<script type=\"text/javascript\" src=\"" << *iter <<"\"> </script>\n";
      htmlStream << "</head>\n";
    }

    void pushInfo(const std::string & text)
    {
      htmlStream << text;
    }

    std::string getDoc()
    {
      return htmlMarkup("html", htmlStream.str());
    }

  private:
    std::ostringstream htmlStream;
  };

  /// Class to draw with the JSXGraph library in HTML page.
  class JSXGraphWrapper
  {
  public:
    JSXGraphWrapper() {
      cpt = 0;
    }

    void reset()
    {
      stream.str("");
      cpt = 0;
    }

    void init(const std::string & sGraphName, int W, int H)
    {
      reset();
      stream
        << "\n"
        << "<div id=\"" << sGraphName << "\" class=\"jxgbox\" style=\"width:"<< W << "px; height:" << H <<"px;\"></div>\n"
        << "<script type=\"text/javascript\">\n"
        << "var board = JXG.JSXGraph.initBoard('"<< sGraphName <<"', {boundingbox: [-10, 10, 10, -10], axis:true,showCopyright:false});\n"
        << "board.suspendUpdate();\n";
    }

    template<typename T>
    void setViewport(const std::pair< std::pair<T,T>, std::pair<T,T> > & range)
    {
      stream
        << "board.setBoundingBox(["
        << range.first.first << ","<< range.second.second <<","
        << range.first.second << ","<< range.second.first <<"]);\n";
    }

    void addLine(double x0, double y0, double x1, double y1, std::string color ="00ff00")
    {
      size_t index0 = cpt++;
      size_t index1 = cpt++;
      stream
        <<"var p"<<index0<<" = board.create('point',["<<x0<<","<<y0<<"], {fixed:true});\n"
        <<"var p"<<index1<<" = board.create('point',["<<x1<<","<<y1<<"], {fixed:true});\n"
        <<"var li = board.create('line',[p"<<index0<<",p"<<index1<<"], "
        <<"{strokeColor:'"<< color <<"',strokeWidth:2});\n";
    }

    template<typename T, typename T2>
    void addXYChart(const std::vector<T> & vec_x, const std::vector<T2> & vec_y,
      std::string stype)
    {
      size_t index0 = cpt++;
      size_t index1 = cpt++;

      stream.precision(5);
      stream.setf(std::ios::fixed,std::ios::floatfield);   // floatfield set to fixed

      stream << "var data"<< index0<<"= [";
      copy(vec_x.begin(), vec_x.end(), std::ostream_iterator<T>(stream, ","));
      stream << "];\n";
      stream << "var data"<< index1<<"= [";
      copy(vec_y.begin(), vec_y.end(), std::ostream_iterator<T>(stream, ","));
      stream << "];\n";
      std::ostringstream osData;
      osData<<"[ data" <<index0<<","<<"data"<<index1<<"]";

      stream << "board.createElement('chart', "
        <<osData.str()
        <<", {chartStyle:'"<< stype <<"',labels:"<<osData.str()<<"});\n";
    }

    template<typename T>
    void addYChart(const std::vector<T> & vec_y, std::string stype)
    {
      size_t index0 = cpt++;

      stream << "var data"<< index0<<"= [";
      copy(vec_y.begin(), vec_y.end(), std::ostream_iterator<T>(stream, ","));
      stream << "];\n";
      stream << "board.createElement('chart', "
        <<"data"<<index0
        <<", {chartStyle:'"<< stype <<"',labels:"<<"data"<<index0<<"});\n";
    }

    void UnsuspendUpdate()
    {
       stream << "board.unsuspendUpdate();\n";
    }
    void close()
    {
      stream << "</script>\n";
    }

    std::string toStr() const
    {
      return stream.str();
    }
  private:
    std::ostringstream stream;
    size_t cpt; //increment for variable
  };
} // namespace htmlDocument

#endif // MIMATTE_HTML_DOC_H
