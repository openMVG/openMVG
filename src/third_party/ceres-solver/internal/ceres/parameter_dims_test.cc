//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: jodebo_beck@gmx.de (Johannes Beck)

#include "ceres/internal/parameter_dims.h"

#include <gtest/gtest.h>

#include <type_traits>
#include <utility>

namespace ceres {
namespace internal {

// Static parameter dims unit test
static_assert(
    std::is_same<StaticParameterDims<4, 2, 1>::Parameters,
                 std::integer_sequence<int, 4, 2, 1>>::value == true,
    "Unit test of type 'parameters' for static parameter dims failed.");

static_assert(StaticParameterDims<4, 2, 1>::kIsValid == true,
              "Unit test of is valid for static parameter dims failed.");
static_assert(StaticParameterDims<4, 2, 1>::kIsDynamic == false,
              "Unit test of is dynamic for static parameter dims failed.");
static_assert(StaticParameterDims<4, 2, 1>::kNumParameterBlocks == 3,
              "Unit test of number of parameter blocks for static parameter "
              "dims failed.");
static_assert(
    StaticParameterDims<4, 2, 1>::kNumParameters == 7,
    "Unit test of number of parameters for static parameter dims failed.");

// Dynamic parameter dims unit test
static_assert(DynamicParameterDims::kIsValid == true,
              "Unit test of is valid for dynamic parameter dims failed.");
static_assert(DynamicParameterDims::kIsDynamic == true,
              "Unit test of is dynamic for dynamic parameter dims failed.");
static_assert(DynamicParameterDims::kNumParameterBlocks == 0,
              "Unit test of number if parameter blocks for dynamic parameter "
              "dims failed.");
static_assert(
    DynamicParameterDims::kNumParameters == 0,
    "Unit test of number of parameters for dynamic parameter dims failed.");

TEST(ParameterDims, GetDims) {
  constexpr int N0 = 3;
  constexpr int N1 = 4;
  constexpr int N2 = 2;

  StaticParameterDims<N0, N1, N2> params;
  EXPECT_EQ(N0, params.GetDim(0));
  EXPECT_EQ(N1, params.GetDim(1));
  EXPECT_EQ(N2, params.GetDim(2));
}

TEST(ParameterDims, GetUnpackedParameters) {
  constexpr int N0 = 3;
  constexpr int N1 = 4;
  constexpr int N2 = 2;

  using ParameterDims = StaticParameterDims<N0, N1, N2>;

  std::array<double, ParameterDims::kNumParameters> packed_parameters{};
  std::array<double*, 3> unpacked_parameters =
      ParameterDims::GetUnpackedParameters(packed_parameters.data());

  EXPECT_EQ(packed_parameters.data(), unpacked_parameters[0]);
  EXPECT_EQ(packed_parameters.data() + N0, unpacked_parameters[1]);
  EXPECT_EQ(packed_parameters.data() + N0 + N1, unpacked_parameters[2]);
}

}  // namespace internal
}  // namespace ceres
