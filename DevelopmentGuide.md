# Developing the VBR calculator

You are welcome to extend the VBR calculator in any way you see fit. This guide is for those who want to do so following the existing framework or for those who want to contribute to the github repo.

## Working with github

To contribute to the github repo:

## Adding a new method
To add a new method:

1. Open the corresponding parameter file in `vbr/4_VBR/VBR_v0p95/params` and then:
  * add the new method name to the `params.possible_methods` cell array
  * add an `elseif` catch for the new method name
  * within the `elseif`, set `param.func_name` to the name of the matlab function that you will write for the new method, e.g.,
  ```
  param.func_name='new_vbr_method'
  ```
  * set any other values that the parameter will need.
2. Create a new file in `vbr/4_VBR/VBR_v0p95/functions` for your new function with the name from `param.func_name`. Using the above example, that would be `new_vbr_method.m`.
3. Write your new method function. The function must have the `VBR` struct as input and output:
```
function [VBR] = new_vbr_method(VBR)
```
The remainder of the function is where you write whatever calculations are appropriate. The VBR structure will come in with all the state variables and parameter value. State variables are accessed with, e.g., ```VBR.in.SV.T``` or ```VBR.in.SV.phi```.
4. To return the results of your function, modify the `VBR.out` structure appropriately, e.g., ```VBR.out.method_type.method_name.result = result;```
where `method_type` is the same as the parameter file that you modified (`anelastic`,`elastic` or `viscous`) and `method_name` is the name you added to `params.possible_methods`

To use your new method, simply add the new method name to the `methods_list`, before you call `VBRspine`, e.g.:
```
VBR.in.method_type.methods_list={'method_name'}
```
where `method_type` is `anelastic`,`elastic` or `viscous` and `method_name` is your new method.
