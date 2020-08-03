# coverage-measure
A measure of correspondence coverage for shape correspondence and registration methods witten for MATLAB. More information about the dataset is available on the [project page](http://robertodyke.com/shrec2020/index2.html).

**Note:** please see `LICENSES.md` for the relevant licensing information about code owned by other individuals.

## Usage
* `compute_coverage.m` this is the main code for computing the coverage measure. Results should be stored in `*.mat` files with the file names `sourcename_targetname.mat` as either `.corr` for vertex-to-vertex correspondences or `.baryc_corr` for vertex-to-barycentric correspondences. Modify the directories specified to point toward the correct ones used.
* `extract_corrs.m` this is for registration methods where the original output is a deformed source model. This code projects the deformed vertices onto the surface of the target shape, saving the output in the necessary format.

If using `toolbox_fast_marching/` to compute geodesics as per default, `toolbox_fast_marching/compile_mex.m` may need to be run first. _(Optional)_ A faster variant of the fast marching method may be found here: https://github.com/orlitany/3D_shapes_tools/. Please refer to the comments in the code of `compute_coverage.m` to easily incorporate this method.

Barycentric coordinates should be stored in a Nx4 matrix. Each row corresponds to a row in the source shape. The first column points to the face that contains the corresponding point, use `-1` to indicate that there is no correspondence for a particular vertex. The subsequent three columns specify the position within a given face with respect to the three vertices, weights must add up to `1`. Example:
```
[ 4  0.2  0.5  0.3
 10  0.7  0.1  0.2
 -1  0.0  0.0  0.0
  7  0.2  0.4  0.4
  3  0.9  0.1  0.0
  5  1.0  0.0  0.0]
```

This work is the product of the following paper:
> Dyke, RM, et al. Shape Correspondence with Non-Isometric Deformations. Computers & Graphics 2020. (accepted)
