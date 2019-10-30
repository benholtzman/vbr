# Developing the VBR calculator

You are welcome to extend the VBR calculator in any way you see fit. This guide is for those who want to do so following the existing framework or for those who want to contribute to the github repository.

## git & github workflow

The github repository follows a stable - master - feature branch framework. The master branch is considered the development branch. Modifications to master branch are made through feature branches. The following workflow outlines the steps for creating and developing a feature branch that eventually merges with master. Only feature branches that are nominally working should be merged into the master branch.

### scope of feature branches

### developing a feature branch

The following assumes that you have already cloned the VBR repository and have a terminal open in the VBR repository's directory `vbr/`

1. **Initial State**: make sure you are on master branch and up to date:
  ```
  git checkout master
  git pull
  ```
2. **New Branch**: create the new feature branch (or branch for fixing a bug):  
  ```
  $git checkout -b new_branch
  ```
where new_branch is the name of your branch. Keep branch names as short and as descriptive as possible.

3. **Push to remote (optional)**: If your branch will take a while to develop, if you want others to be able to checkout or contribute to your branch, or if you want to use multiple computers to develop your branch, push your branch to the remote repository with `--set-upstream` :
  ```
  git push --set-upstream origin new_branch
  ```
This command sets the remote branch that your local `new_branch` will track. Any `git push` or `pull` will now automatically sync with the remote `new_branch` on github.

4. **Develop your branch**: develop as normal on the new branch, adding commits as you see fit.

5. **Test your branch**: If your feature branch adds new functionality to the `vbr` directory, you should occasionally run the test functions in `vbr/testing` (see the README there, `vbr/testing/README.md`) during development and add new test functions for your new features. If your new feature branch is a self contained project in `Projects`, new test functions are not required. In either case, before merging back to master, please run the full test.

### merging your new feature branch back into master

### merging latest master branch into your feature branch
If your feature branch needs to use updates from other feature branches that have been merged into the master branch, then you can pull those changes into your new branch in several ways.

### conflicts

### style guide & helpful git tips:

In case you're new to git or developing VBR, here are some helpful tips!

**commits**: When committing code, please keep the as a brief description. When a longer description is useful, keep the first line of the commit as a brief description and add a longer description on the third line (to do this, it's easier to use `git commit` and edit the commit in your text editor rather than a single line commit, `git commit -m "commit message"`).

**stashing**: If you need to swtich branches but aren't ready to commit changes to your code, you can stash your uncommitted changes. `git stash` will store uncommitted changes on current branch, `git stash apply` will restore those uncommited changes on your current branch.

**hot fixes**: for very quick bug fixes (typos, etc.), editing on the master branch is OK. If you're not sure whether you need a new branch, you probably need a new branch. 

## Adding a new VBR core method
To add a new method to the VBR core:

1. Open the corresponding parameter file in `vbr/vbrCore/params` and then:
  * add the new method name to the `params.possible_methods` cell array
  * add an `elseif` catch for the new method name
  * within the `elseif`, set `param.func_name` to the name of the matlab function that you will write for the new method, e.g.,
  ```
  param.func_name='new_vbr_method'
  ```
  * set any other values that the parameter will need.
2. Create a new file in `vbr/vbrCore/functions` for your new function with the name from `param.func_name`. Using the above example, that would be `new_vbr_method.m`.
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
