



![Alt text](images/title.png)

Contourklip is a self-contained single-header C++ library for boolean operations on multipolygons/contours ("polygon clipping"), where contours (paths) can consist of line segments *and* cubic bezier curves. 

The library focuses on correctness and a straightforward API. It is a modification and extension of the algorithm from the paper by Martínez[^1]. It supports the 5 common boolean operations, namely *union, intersection, difference, xor (symmetric difference), divide*. Please note that this library is in a relatively early state.



------



## Installation

Contourklip needs no dependencies and consists of a single header, simply add the header in `single_include` to your project. 

## Usage

This section is currently considered to be the documentation. Here is a basic example to get started, which corresponds to the computation shown in the title picture:

```c++
#include <iostream>
#include "contourklip.hpp"

int main() {
    contourklip::Contour contour1{{0, 100}};
    contour1.push_back({50, 100});
    contour1.push_back({77.5, 100}, {100, 77.5}, {100, 50});
    contour1.push_back({100, 22.5}, {77.5, 0}, {50, 0});
    contour1.push_back({0, 0});
    contour1.close();

    contourklip::Contour contour2{{150, 25}};
    contour2.push_back({100, 25});
    contour2.push_back({72.3, 25}, {50, 47.3}, {50, 75});
    contour2.push_back({50, 102.5}, {72.3, 125}, {100, 125});
    contour2.push_back({150, 125});
    contour2.close();

    std::vector<contourklip::Contour> shape1{contour1};
    std::vector<contourklip::Contour> shape2{contour2};
    std::vector<contourklip::Contour> result{};

    if(contourklip::clip(shape1, shape2, result, contourklip::INTERSECTION)){
        std::cout << "clipping operation succeeded\n";
    }

    for (const auto &contour: result) {
        std::cout << "contour:\n";
        std::cout << contour;
    }
    return 0;
}
```

Output:

```
clipping operation succeeded
contour:
(50, 75)
(50, 49.5052) (68.8907, 28.5849) (93.4955, 25.4155)
(97.6328, 32.6861) (100, 41.08) (100, 50)
(100, 75.312) (80.9379, 96.388) (56.4602, 99.5815)
(52.3456, 92.3116) (50, 83.9187) (50, 75)
```

In this case there is 1 contour in the result, and one can see that it indeed makes sense compared to the expected result (namely a "leaf" shape).


Let's take a closer look:

- Here, a multipolygon is indeed simply a vector of `Contour`s.

- To compute a boolean operation we just have to call the following function, which looks quite intuitive:

  ```c++
  bool clip(const std::vector<Contour> &a, const std::vector<Contour> &b,
                    std::vector<Contour> &result, BooleanOpType clippingop)
  ```

- The `BooleanOpType` should be one of `UNION, INTERSECTION, DIFFERENCE, XOR, DIVIDE`, as one would expect.

- The function will store the result in `out`, and returns if the clipping operation succeeded. The function will only append to the result, and not modify or remove contours which are already in the result vector.

- If one has mutipolygons which contain only 1 contour anyway, there's a corresponding overload for convenience:

  ```c++
  bool clip(const Contour &a, const Contour &b,
                    std::vector<Contour> &result, BooleanOpType clippingop)
  ```

  Note that the result remains a vector of contours. Consider for example the union of two nonintersecting polygons. 

IMPORTANT: at this point it is worth mentioning that if you are getting unexpected output, make sure to read the "assumptions and limitations" section carefully. 

##### Contours & ContourComponents

A `Contour` is a random-access container which represents a contour, by containing `ContourComponent`s. 

A `ContourComponent` represents a line segment which is connected to the previous point (1 point is needed), or a cubic bezier which is connected to the previous point (3 points are needed). So, the `ContourComponent` has 3 points and whether the `ContourComponent` is one or the other is determined by an `enum ComponentType {LINE = 0, CUBIC_BEZIER = 1}` . (Indeed there's no subclassing for line and bezier). If it is a line, the other 2 points are not meaningful.

 The following methods are provided:

- `explicit ContourComponent(const Point2d &pLast)` initializes it so that it represents a line
- `ContourComponent(const Point2d &c_1, const Point2d &c_2, const Point2d &p)` initializes it so that it represents a cubic bezier.
- `ComponentType segment_shape() const` returns the enum associated with this instance (`LINE` or `CUBIC_BEZIER`)
- `Point2d point() const `, `Point2d &point()` returns endpoint/ last point of the implied segment when connected to a previous point (in a contour).
- `Point2d c1() const`, `Point2d& c1()` returns the first bezier control point. Meaningless if `this->segment_shape() == LINE`. 
- `Point2d c2() const`, `Point2d& c2()` returns the second bezier control point. Meaningless if `this->segment_shape() == LINE`. 

While `ContourComponent`s need not be used to initialize contours, they are useful to work with the result. 

Then, `Contour` provides the following methods:

- `explicit Contour(const Point2d &start)` Constructor with starting point. This is what should be used by default. 
  
- `Contour(const Point2d &p0, const Point2d &c1)` initializes a contour from a line segment, for convenience.

- `Contour(const Point2d &p0, const Point2d &c1, const Point2d &c2, const Point2d &p3) ` Initializes a contour from a cubic bezier, for convenience.

-  `Contour()` trivial constructor. Note that this can lead to degenerate contours if then the first component pushed back represents a cubic bezier, because then the starting point is missing. Therefore the first construct with a point is preferrable. 

- `void push_back(const Point2d &p)` appends a point, which creates a line segment with the previous point. 

- `void push_back(const Point2d &c2, const Point2d &p3, const Point2d &p)` appends a cubic bezier which is connected to the previous point.

- `void push_back(const ContourComponent &start)` appends a `ContourComponent `. It is preferable to use the other methods. 

- `ContourComponent operator[](const size_t idx) const` , `ContourComponent& operator[](const size_t idx)` : returns the component at index idx. 

- `std::size_t size() const` returns the number of `ContourComponent`s. 

- `Point2d front_point() const` returns the starting point.

- `Point2d back_point() const` returns the last segment point.

- `ContourComponent &front()`, `ContourComponent front() const`: returns the first `ContourComponent`

- `ContourComponent &back()`, `ContourComponent back() const`: returns the last `ContourComponent`

- `bool is_closed() const` Returns true iff contour is closed, that is, the last segment point corresponds to the starting point.

- `void close()` Closes the contour if not closed.

- `void reverse()` Reverses the contour direction in-place.

- `template<ComponentType T, typename Consumer>  void forward_segments(Consumer &out) const` pass all line segments or all cubic beziers (according to `ComponentType`) to a consumer callback `Consumer`, by passing 2 points for a line, or 4 point for a cubic bezier. The following will print all cubic beziers of a contour:

  ```c++
  contourklip::Contour c1{{0, 0}};
  c1.push_back({1, 0});
  c1.push_back({2, 0}, {2, 1}, {1, 1});
  c1.push_back({0, 1});
  c1.close();
  auto print_lines = [](const Point2d& a, const Point2d& b){
    std::cout << "line: " << a << " " << b << "\n";
  };
  auto print_curves = [](const Point2d& p0, const Point2d& c1, const Point2d& c2, const Point2d& p3){
    std::cout << "bezier: " << p0 << " " << c1 << " " << c2 << " " << p3 << "\n";
  };
  c1.forward_segments<contourklip::LINE>(print_lines);
  c1.forward_segments<contourklip::CUBIC_BEZIER>(print_curves);
  ```

  Output:

  ```
  line: (0, 0) (1, 0)
  line: (1, 1) (0, 1)
  line: (0, 1) (0, 0)
  bezier: (1, 0) (2, 0) (2, 1) (1, 1)
  ```

  

-  `auto begin() const`, `auto begin()`, `auto end()`,  `const auto end()` : iterator pairs, iterating over the `ContourComponent`s. 

Summing up, the following can be used to actually retrieve the points:

```c++
for (const auto &contour: result) {
  for (const auto &seg: contour) {
    switch (seg.segment_shape()) {
      case contourklip::LINE:
        //process line, extract point coordinates
        std::cout << "(" << seg.point().x() << ", " <<seg.point().y() << ")";
        break;
      case contourklip::CUBIC_BEZIER:
        //process curve
        std::cout << seg.c1()
          <<", " << seg.c2() << ", " << seg.point();
        break;
    }
    std::cout << "\n";
  }
}
```

##### Slightly advanced usage

As pointed out in the introduction, more fine-grained usage is quite limited. In particular, `double` is the default numeric type, and the given classes need to be used, for example there's no `point` template. 

The function discussed previously is just a wrapper around the ` contourklip::Polyclip` class: 

```c++
 template<typename Orient2dFunc = LeftOfLine, typename CollinearFunc = IsCollinear>
    class PolyClip {
    PolyClip(const std::vector<Contour> &a, const std::vector<Contour> &b,
             std::vector<Contour> &result, BooleanOpType clippingop, Config c = {});
		bool success();
		void compute();
}
```

The workflow is thus as follows: we initialize an instance with appropriate parameters, call `compute()` and then check `success()`. 

There are 2 template parameters for geometry predicates: `Orient2dFunc`, to check if a point is on the left of a directed segment defined by 2 other points, and `CollinearFunc`, to check if 3 points are collinear. These default to an implementation which might not be robust with respect to roundoff errors in unlikely degenerate cases. 

Additionally, an optional `contourklip::Config` struct can be passed. It has the following fields:

- `postprocess`: whether to perform any post-processing on a result contour. Concretely this means to split it into simpler contours if possible (i.e. if there are overlapping points) , and to remove points which are not needed (see next field). Defaults to `true`. 
- `remove_collinear`: wether to remove successive collinear points in the result contours. Defaults to `true`. Implied to be false when post-process is false. 
- `fail_on_approx_equal`: whether to set success to false if approximately equal points are detected. This is motivated by the fact that it likely indicates a numerical error. Defaults to `true`.
- `approx_equal_tol`: the absolute tolerance to decide if 2 points are approximately equal. Defaults to `1e-8`. 

We illustrate this with the following example, which in this case ultimately performs the exact same operation as in the first example. Note that one needs c++>= 20 to allow default-constructible non-capturing lambdas. Otherwise one needs a functor, i.e. a struct with a `bool operator()` overload. 

```c++
#include <iostream>
#include "contourklip.hpp"

int main() {
  
    // ...contour initialization same as before
  
    std::vector<contourklip::Contour> shape1{contour1};
    std::vector<contourklip::Contour> shape2{contour2};
    std::vector<contourklip::Contour> result{};

    //callback for determining if a is on the left of the segment p0, p1.
    // here it just calls the default callback.
    auto above =
            [](const contourklip::Point2d &p0,
                    const contourklip::Point2d &p1, const contourklip::Point2d &a) -> bool{
        return contourklip::detail::LeftOfLine{}(p0, p1, a);
    };

    // callback for determining if 3 points are collinear.
    // Again this just calls the default implementation.
    auto collinear =
            [](const contourklip::Point2d &p0,
                    const contourklip::Point2d &p1, const contourklip::Point2d &p2) -> bool{
        return contourklip::detail::IsCollinear{}(p0, p1, p2);
    };

    contourklip::Config c;
    c.postprocess = false;

    contourklip::PolyClip<decltype(above), decltype(collinear)> clip{shape1, shape2, result,
                                                                     contourklip::INTERSECTION, c};
    clip.compute();

    if(clip.success()){
        std::cout << "clipping operation succeeded\n";
    }
    return 0;
}
```



## Assumptions and Limitations

##### Definitions

- A segment is either a line segment or a cubic bezier, and in the second case we also say bezier segment. 

- A contour is an ordered list of connected segments. 

- A multipolygon is a set of contours which defines a shape. While it makes sense to think of it as an svg path, it is not quite accurate since an svg path has more features. Also note that the term "polygon" is not ideal since that implies there are only lines, but "multipolygon" is more common than "multicontour".

##### Input contour treatment

The algorithm always uses the "evenodd" fill rule when interpreting a multipolygon as a shape (see [here](https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/fill-rule)). One way to think about this is that all closed shapes that appear in the (same) multipolygon are XORed together. Hence, the contour direction (clockwise/counterclockwise) of any input contour is *not* taken into account. 

The algorithm also assumes each contour is closed, since clipping open contours is not well defined. 

##### Shape structure of the output 

Although the visual output of a clipping operation is well-defined, sometimes are multiple ways the output shape can be represented with contours. In general, no assumption is made about the output contour representation, in particular with respect to contour direction. However, an attempt is made to "disentangle" the output by post-processing it. Also note that for two visually identical inputs represented differently in terms of contours, the output representation may not be the same. 

A particular focus was put on the property that the output has clean curve geometry. This means that the output only has as much beziers as needed. 

##### Preconditions

The library makes certain assumptions about the input contours. The following may lead to an unsuccessful computation:

- Two segments (line or bezier) of the *same* contour overlap in more than just one point (i.e. they share a subsegment)
- A line segment has the same startpoint and endpoint
- A bezier segment intersects itself (this might eventually be supported)
- A bezier segment has the same startpoint and endpoint (this is a special case of the previous point)
- two bezier segments share a point of tangency (i.e. they have a common point but without intersecting) which does not corrrespond to a start/endpoint.
- A startpoint or endpoint (of a curve or line segment) is *on* another bezier curve B, unless it's on the startpoint or endpoint of the curve B.
- Input data contains NaN or Inf values. 

Otherwise there are no restrictions, and in particular a contour can intersect itself. 

However, while the algorithm then *should* compute the correct result, it is best to still check wether it succeeded by checking `success()`, in particular because of numerical errors. 

##### Circles, quadratic beziers and other curves

The only curve type that is supported is the cubic bezier curve. As such, not all svg paths can be implemented in a `contourklip::Contour`. However note that quadratic beziers can be trivially elevated to cubic beziers, and circles or ellipses are often approximated with cubic beziers (indeed the approximation is very accurate). 

##### A note about divide

The *divide* operation can be thought of as retrieving all closed shapes that appear when looking at all overlapping contour segments. Here it is defined as the concatenation of *xor* and *intersection*. Note that this is not quite equivalent.

##### Numerical considerations

- Clearly, the library actually computes an accurate numerical approximation of a clipping operation, not least because the intersection of beziers has no closed-form solution. After all the purpose of this library is to have an explicit shape representation of the clipping operation result. 
- The implementation might not be robust, but numerical issues leading to an incorrect state should be detected in which case `success()` returns false. Also, since contourklip is self-contained, it does not use a robust geometry kernel such as CGAL. However, as shown in the "advanced usage" section, some predicate functions from eg. CGAL can be used. 

## About the algorithm

Any clipping algorithm has to do the following:

- calculate all segment intersections
- determine which (sub)segments belong to the result
- construct the result

How to do it properly is nontrivial. The library uses the code/main idea from the paper by Martínez[^1], whose implementation is available on his website[^2]. Note that the paper is not concerned with beziers but "only" standard line polygons. The code is in the public domain and a few parts have been refactored into the library (with additional permission from the author). Other than that there are some significant differences, for example the original algorithm also computes the intersections during the plane sweep, which leads to very complicated (but more efficient) code. Instead we compute all pairwise intersections. Last but not least, support for cubic beziers has been added. 

To compute intersections of cubic bezier curves and generally do bezier-related operations, techniques from here[^3] have been used.  

Those working with very large *line* polygon data definitely may want to look at the original paper code instead (or use another library of course, such as Boost.Geometry). 



[^1]: https://doi.org/10.1016/j.advengsoft.2013.04.004
[^2]: https://www4.ujaen.es/~fmartin/bop12.zip
[^3]:  https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=1000&context=facpub
