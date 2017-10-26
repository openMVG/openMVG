// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014-2016 Yongchao Xu, Pascal Monasse, Thierry Géraud, Laurent Najman
// Copyright (c) 2016 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/tbmr/tbmr.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_container.hpp"

#include <numeric>

namespace openMVG
{
namespace features
{
namespace tbmr
{
  namespace internal
  {
    static
    inline
    unsigned int
    zfindroot(image::Image<unsigned int>& parent, unsigned int p)
    {
      if (parent[p] == p)
        return p;
      else
        return parent[p] = zfindroot(parent, parent[p]);
    }
  }

  template <typename I, typename BinaryFunction = std::less<I>>
  std::vector<unsigned int>
  pixel_indexes_ordering
  (
    const image::Image<I>& input,
    BinaryFunction cmp = BinaryFunction ()
  )
  {
    const unsigned int pixel_count = input.Width()*input.Height();
    std::vector<unsigned int> v(pixel_count);
    std::iota(v.begin(), v.end(), 0);
    std::sort(v.begin(), v.begin() + pixel_count,
      [&input, cmp](unsigned int x, unsigned int y)
      { return cmp(input[x], input[y]); }
    );
    return v;
  }

  std::vector<Vec2i>
  wrt_delta_index
  (
  )
  {
    return { { 0,-1 }, {0,1}, {-1, 0}, {1,0} };
  }

  //for incremental computation of region information
  struct attribute{
    unsigned int area = 0;
    double sum_x = 0;
    double sum_y = 0;
    double sum_xy = 0;
    double sum_xx = 0;
    double sum_yy = 0;

    attribute& operator += (const attribute & rhs)
    {
      area   += rhs.area;
      sum_x  += rhs.sum_x;
      sum_xx += rhs.sum_xx;
      sum_xy += rhs.sum_xy;
      sum_y  += rhs.sum_y;
      sum_yy += rhs.sum_yy;
      return *this;
    }

    void init
    (
      unsigned int area,
      double x,
      double y
    )
    {
      this->area = area;
      sum_x = x;
      sum_y = y;
      sum_xy = x * y;
      sum_xx = x * x;
      sum_yy = y * y;
     }
  };

  // Template instantiation for WHITE features
  template
  void Extract_tbmr
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::AffinePointFeature> & features,
    std::less<unsigned char> cmp,
    const unsigned int minimumSize,
    const double maximumRelativeSize
  );

  // Template instantiation for DARK features
  template
  void Extract_tbmr
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::AffinePointFeature> & features,
    std::greater<unsigned char> cmp,
    const unsigned int minimumSize,
    const double maximumRelativeSize
  );

  template <typename Ordering>
  void Extract_tbmr
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::AffinePointFeature> & features,
    Ordering cmp,
    const unsigned int minimumSize,
    const double maximumRelativeSize
  )
  {
    //construct the Max/Min-tree based on union-find [Berger 2007 ICIP] + rank
    //and compute incrementally the moments information during tree construction
    //--------------------------------------------------------------------------

    //sort the pixels in the order defined by cmp
    const std::vector<unsigned int> S = pixel_indexes_ordering(ima, cmp);

    std::vector<unsigned int> parent(ima.Width() * ima.Height());
    std::vector<unsigned int> root(ima.Width() * ima.Height());
    std::vector<unsigned int> rank(ima.Width() * ima.Height(), 0);
    std::vector<bool> dejaVu(ima.Width() * ima.Height(), false);
    std::vector<attribute> imaAttribute(ima.Width() * ima.Height());
    image::Image<unsigned int> zpar(ima.Width(), ima.Height());

    const std::vector<Vec2i> offsets(wrt_delta_index());

    for (int i = S.size()-1; i >= 0; --i)
    {
      const unsigned int p = S[i];
      // make set
      {
        parent[p] = p;
        zpar[p] = p;
        root[p] = p;
        dejaVu[p] = true;
        imaAttribute[p].init(1, p % ima.Width(), p/ima.Width());
      }

      unsigned int x = p; // zpar of p
      for (unsigned int k = 0; k < offsets.size(); ++k)
      {
        const Vec2i point_p (p % ima.Width(), p / ima.Width());
        const Vec2i point_q = point_p + offsets[k];
        const unsigned int q = point_q[1] * ima.Width() + point_q[0];
        if (ima.Contains(point_q[1],point_q[0]) && dejaVu[q])
        {
          const unsigned int r = internal::zfindroot(zpar, q);
          if (r != x)
          { // make union
            parent[root[r]] = p;
            //accumulate information
            imaAttribute[p] += imaAttribute[root[r]];

            if (rank[x] < rank[r])
            {
              //we merge p to r
              zpar[x] = r;
              root[r] = p;
              x = r;
            } else
            if (rank[r] < rank[p])
            {
              // merge r to p
              zpar[r] = p;
            } else
            {
              // same height
              zpar[r] = p;
              rank[p] += 1;
            }
          }
        }
      }
    }
    //end of Min/Max-tree construction
    //--------------------------------------------------------------------------

    // canonization
    for (unsigned int p: S)
    {
      const unsigned int q = parent[p];
      if (ima[parent[q]] == ima[q])
        parent[p] = parent[q];
    }

    //TBMRs extraction

    /* small variant of the given algorithm in the paper. For each
    critical node having more than one child, we check if the
    largest region containing this node without any change of
    topology is above its parent, if not, discard this critical
    node */

    /* note also that we do not select the critical nodes themselves
     * as final TBMRs */

    //--------------------------------------------------------------------------
    std::vector<unsigned int> numSons(ima.Width() * ima.Height(), 0);
    std::vector<unsigned int> vecNodes(imaAttribute[S[0]].area);
    unsigned int numNodes = 0;

    //leaf to root propagation to select the canonized nodes
    for (int i = S.size()-1; i >= 0; --i)
    {
      const unsigned int p = S[i];
      if (parent[p] == p || ima[p] != ima[parent[p]]) {
        vecNodes[numNodes++] = p;
        if (imaAttribute[p].area >= minimumSize)
          ++numSons[parent[p]];
      }
    }

    std::vector<bool> isSeen(ima.Width() * ima.Height(), false);

    //parent of critical leaf node
    std::vector<bool> isParentofLeaf(ima.Width() * ima.Height(), false);
    for (int i = 0; i < vecNodes.size(); ++i)
    {
      const unsigned int p = vecNodes[i];
      if (numSons[p] == 0 && numSons[parent[p]] == 1)
        isParentofLeaf[parent[p]] = true;
    }

    const unsigned int maxArea = maximumRelativeSize * S.size();
    unsigned int numTbmrs = 0;

    std::vector<unsigned int> vecTbmrs(numNodes);
    for (const unsigned p : vecNodes)
    {
      if (numSons[p] == 1 && !isSeen[p] && imaAttribute[p].area <= maxArea)
      {
        unsigned int num_ancestors = 0;
        unsigned int pt = p;
        unsigned int po = pt;
        while (numSons[pt] == 1 && imaAttribute[pt].area <= maxArea)
        {
          isSeen[pt] = true;
          ++num_ancestors;
          po = pt;
          pt = parent[pt];
        }
        if (!isParentofLeaf[p] || num_ancestors > 1)
        {
          vecTbmrs[numTbmrs++] = po;
        }
      }
    }
    //end of TBMRs extraction
    //------------------------------------------------------------------------

    //compute best fitting ellipses
    //------------------------------------------------------------------------
    for (unsigned int i = 0; i < numTbmrs; ++i)
    {
      const unsigned int p = vecTbmrs[i];
      const double x = imaAttribute[p].sum_x / (double)imaAttribute[p].area;
      const double y = imaAttribute[p].sum_y / (double)imaAttribute[p].area;

      const double i20 = imaAttribute[p].sum_xx - imaAttribute[p].area*x*x;
      const double i02 = imaAttribute[p].sum_yy - imaAttribute[p].area*y*y;
      const double i11 = imaAttribute[p].sum_xy - imaAttribute[p].area*x*y;
      const double n = i20*i02 - i11*i11;

      if (n != 0)
      {
        const double a = i02/n * (imaAttribute[p].area-1)/4;
        const double b = -i11/n * (imaAttribute[p].area-1)/4;
        const double c = i20/n * (imaAttribute[p].area-1)/4;

        const features::AffinePointFeature affineFP (x, y, a, b, c);

        // Check feature validity  (avoid tiny and thick ellipses)
        const double lMin = std::min(affineFP.l1(), affineFP.l2());
        if (lMin < 1.5)
          continue;
        // TODO: delete the one that collide with the border !!
        features.emplace_back(affineFP);
      }
    }
  }

} // namespace tbmr
} // namespace features
} // namespace openMVG
