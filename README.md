# MAC
### Multidimensional Array Containers.

This is a tiny modern Fortran library to create and manipulate arrays of any rank. The library is meant to serve as a building block for codes that expand on "rank-agnostic" programming. See, e.g., [this](https://fortran-lang.discourse.group/t/examples-for-the-new-rank-agnostic-array-notation/1376/8) thread.

# API

The following derived types are defined:
``` fortran
type, public :: container_specifier
type, extends(container_specifier), public :: container
```
The former is an array handle, which stores the rank, shape, bounds and the [layout](https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays) of the multidimensional array. The latter is an extension for the handle, which provides the storage for the actual array.

## `type(container_specifier) :: a`
Procedure | Description | Parameters
--- | --- | ---
Constructor subroutine. <br /> Usage: <br /> `call a%specify(dimension_specifier, [lower_bounds, layout])` | Initializes the array handle. | `integer, intent(in) :: dimension_specifier(:)` : 1D array containing the shape of the multidimensional array. <br />`integer, intent(in) :: lower_bounds(:)` (optional) : 1D array containing the shape of the multidimensional array. Default is an array of 1-s. <br />`character(len = *), intent(in) :: layout` (optional) : specifies which dimension should contain contiguous indexes upon traversal of the array. The options are `"F"`: leftmost dimension; and `"C"`: rightmost dimension. Default is `"F"`. [See also](https://fortran-lang.org/learn/best_practices/multidim_arrays/). <br /> `type(container_specifier), intent(out) :: a`: initialized array handle.
Size handle. <br /> Usage: <br /> `sz = a%size()`| Returns the size of the represented array.| `type(container_specifier), intent(in) :: a`: initialized array handle. <br /> `integer, intent(out) :: sz`: size of the represented array.
Indexing functions. <br /> Usage: <br /> `[arr], {mem} = a%ind([mem], {arr})` | Shifts between memory layout representation and array layout representation. | `type(container_specifier), intent(in) :: a`: initialized array handle. <br /> Given `mem`: <br /> `integer, intent(in) :: mem`: integer referencing a memory layout adress. Must be in the range `mem`$\in[1,$ `a%size()`$]$. <br /> `integer, intent(out) :: arr(:)`: 1D integer array referencing the corresponding array layout adress. The number of components is equal to the rank of the represented array. Each component `arr(i)` is in the range `arr(i)`$\in[$`lb(i)`$,$ `lb(i) + spc(i) - 1`$]$. Where `lb(i)` and `spc(i)` are the lower bounds and number of elements for dimension `i` of the represented array, respectively. <br /> Given `arr(:)`: <br /> `integer, intent(in) :: arr(:)`: 1D integer array referencing an array layout adress. Each component `arr(i)` must be the range `arr(i)`$\in[$`lb(i)`$,$ `lb(i) + spc(i) - 1`$]$. <br /> `integer, intent(out) :: mem`: integer array referencing the corresponding memory layout adress. Will be in the range `mem`$\in[1,$ `a%size()`$]$.
Property accessor. <br /> Usage: <br /> `call a%get(property, val)` | Retrieves the value of the requested property. | `type(container_specifier), intent(in) :: a`: initialized array handle. <br /> `character(len=*), intent(in) :: property`: requested property, options are `"size"`, `"layout"`, `"dimension"`, `"shape"`, `"lower_bounds"`, `"upper_bounds"`. <br /> `integer, intent(out) :: val[(:)]`: value of the requested property. Returns: the size, the data layout (`val = 0` is `"F"` layout and `val = 1` is `"C"` layout), the dimension or rank, the shape, and the lower and upper bounds of the represented array, respectively.

## `type(container) :: b`

Procedure(*) | Description | Parameters
--- | --- | ---
Constructor subroutine. <br /> Usage: <br /> `call b%construct(container_type, dimension_specifier, [lower_bounds, layout])` | Initializes the array handle and array. | `character(len=*), intent(in) :: container_type`: data type to allocate storage for. Options are `logical`, `integer`, `real`, `complex`, `real_dp`, `complex_dp` where the last two are double precision kinds. <br />`integer, intent(in) :: dimension_specifier(:)` : 1D array containing the shape of the multidimensional array. <br />`integer, intent(in) :: lower_bounds(:)` (optional) : 1D array containing the shape of the multidimensional array. Default is an array of 1-s. <br />`character(len = *), intent(in) :: layout` (optional) : specifies which dimension should contain contiguous indexes upon traversal of the array. The options are `"F"`: leftmost dimension; and `"C"`: rightmost dimension. Default is `"F"`. <br /> `type(container), intent(out) :: b`: initialized array.
Setter subroutine. <br /> Usage: <br /> `call b%set(value [, at = [mem], {arr}])`. | Sets the value of the represented array (component). | `type(container), intent(in) :: b`: initialized array. <br /> `$kind$, intent(in) :: value`: variable to set the data from. The type and kind of this value must be compatible with the container's type and kind. <br /> `integer, intent(in) :: mem`<br /> or <br /> `integer, intent(in) :: arr(:)`<br /> The memory or array layout index where the value is to be set. The restrictions on `mem` and `arr(:)` are the same as for the indexing function. If `at` is not present it sets the value of all components of the represented array.
Getter subroutine. <br /> Usage: <br /> `call b%get(value [, at = [mem], {arr}])`. | Retrieves the value of the represented array component. | `type(container), intent(in) :: b`: initialized array. <br /> `$kind$, intent(out) :: value`: variable to retrieve the data to. The type and kind of this value must be compatible with the container's type and kind. <br /> `integer, intent(in) :: mem`<br /> or <br /> `integer, intent(in) :: arr(:)`<br /> The memory or array layout index of the value is to be retrieved. The restrictions on `mem` and `arr(:)` are the same as for the indexing function.

(*)Additionally to the procedures that apply to `type(container_specifier) :: a`.

Component | Description
--- | ---
Logical storage. <br /> `b%l_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `logical` type.
Integer storage. <br /> `b%i_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `integer` type.
Real storage. <br /> `b%r_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `real` type.
Complex storage. <br /> `b%c_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `complex` type.
Real double precision storage. <br /> `b%rdp_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `real_dp` type.
Complex double precision storage. <br /> `b%cdp_storage(:)` | Memory layout array representing `b`. Will **only** be allocated if the container type has been defined of `complex_dp` type.

# Limitations

The only serious limitation of this approach to rank-agnostic programming is integer overflow. The default integer type is 4-byte, which means that the largest array possible to be represented can have only up to $2^{31} - 1 = 2147483647$ components (or $\approx$ 8 GB for a real array). The library will set a runtime error if this number of components is excedeed upon initailization of the array or array handle by the constructor routines.

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use MAC in your projects. You can add MAC to your project dependencies by including

```
[dependencies]
MAC =
```
to the `fpm.toml` file.
