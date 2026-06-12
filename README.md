# sLORETA Plotting Tools

This folder contains MATLAB tools for running a simple subject-level sLORETA
analysis and saving the results as:

- a five-view brain image
- a ranked brain-region CSV file
- a `.mat` file that can be used later for group-level analysis

The main function is:

```matlab
doPlotsLoreta
```

For condition contrasts, the main function is:

```matlab
doContrastLoreta
```

Most users should start with one of these examples:

```matlab
exampleUseage
exampleContrast
exampleGroup
```

Open `exampleUseage.m` in MATLAB, edit the condition number or time point if
needed, and press Run.

## Quick Start

In MATLAB:

```matlab
cd('/path/to/sLoretta')
exampleUseage
```

The example runs subject 1, condition 1 from a small example data file.

It creates files in the `outputs` folder:

- `rewp_sub01_cond01_sixview.png`
  - the five-view BrainNet ICBM152 brain image
- `rewp_sub01_cond01_slices.png`
  - sagittal, coronal, and axial internal source slices
- `rewp_sub01_cond01_mri_slices.png`
  - sagittal, coronal, and axial source overlays on a standard MRI
- `rewp_sub01_cond01_regions.csv`
  - brain regions ranked from most active to least active
- `rewp_sub01_cond01_source.mat`
  - source-localized data and metadata for later analysis

## What Data Shape Is Expected?

The full group data shape is:

```text
channels x time x conditions x subjects
```

For the beginner examples, one subject has already been extracted into:

```text
exampleSubjectData.mat
```

Inside that file, the variable is:

```text
exampleData
```

Its size is:

```text
63 channels x 400 time points x 2 conditions
```

Data loading is meant to happen outside the main subject-level functions.
This means your variable can be called anything you want. Load it in your
script, then pass the numeric matrix into `doPlotsLoreta` or
`doContrastLoreta`.

## The Main Function

The example function call looks like this:

```matlab
outputs = doPlotsLoreta('exampleSubjectData.mat', ...
    'DataVariable', 'exampleData', ...
    'SubjectIdx', 1, ...
    'ConditionIdx', 1, ...
    'OutputDir', 'outputs', ...
    'OutputPrefix', 'rewp_sub01_cond01');
```

### What Each Input Means

`'exampleSubjectData.mat'`

The `.mat` file containing the EEG/ERP data.

`'DataVariable', 'exampleData'`

The name of the variable inside the `.mat` file.

`'SubjectIdx', 1`

Analyze the first subject.

`'ConditionIdx', 1`

Analyze the first condition.

`'OutputDir', 'outputs'`

Save results in a folder called `outputs`.

`'OutputPrefix', 'rewp_sub01_cond01'`

Start each output filename with `rewp_sub01_cond01`.

The beginner example instead loads a numeric matrix first:

```matlab
load('exampleSubjectData.mat', 'exampleData');

outputs = doPlotsLoreta(exampleData, ...
    'ConditionIdx', 1, ...
    'OutputDir', 'outputs', ...
    'OutputPrefix', 'rewp_sub01_cond01');
```

Both styles are supported. Passing a numeric matrix is the cleaner approach
for a general toolbox because the function does not need to know your file
or variable names.

## Single Condition Or Two Conditions?

Use `doPlotsLoreta` for a single condition.

The input can be either:

```text
channels x time
```

or:

```text
channels x time x conditions
```

If your data includes multiple conditions, set:

```matlab
'ConditionIdx', 1
```

Use `doContrastLoreta` when you want to compare two conditions. The input
must include a condition dimension:

```text
channels x time x conditions
```

The contrast is always:

```text
condition 1 - condition 2
```

The code does not care whether those conditions are wins/losses, targets/
standards, P300 conditions, or anything else. Arrange your matrix so the
condition you want as the positive side of the contrast is condition 1.

## Choosing A Time Point

By default, `doPlotsLoreta` automatically chooses the time point with the
largest average source activity.

For group-level analysis, it is usually better to use the same time point for
every subject.

To force a specific time point, add `TimeIndex`:

```matlab
outputs = doPlotsLoreta('exampleSubjectData.mat', ...
    'DataVariable', 'exampleData', ...
    'SubjectIdx', 1, ...
    'ConditionIdx', 1, ...
    'TimeIndex', 317, ...
    'OutputDir', 'outputs', ...
    'OutputPrefix', 'rewp_sub01_cond01_t317');
```

## Output Files

### Brain Image

The image file ends with:

```text
_sixview.png
```

It shows:

```text
left   top   right
front  back
```

The filename keeps `_sixview` for compatibility with earlier versions, but
the current default figure uses five views and omits the bottom view.

By default the surface is the included BrainNet ICBM152 cortex:

```matlab
'BrainTemplate', 'brainnet'
```

Other built-in choices are:

```matlab
'BrainTemplate', 'brainnet_smoothed'
'BrainTemplate', 'brainstorm'
```

Advanced users can provide a custom BrainNet `.nv` file or MATLAB surface
file:

```matlab
'SurfaceFile', '/path/to/my_surface.nv'
```

MATLAB surface files can contain `Vertices`/`Faces` or FieldTrip-style
`mesh.pos`/`mesh.tri`.

### Internal Slice Image

The slice image ends with:

```text
_slices.png
```

It shows sagittal, coronal, and axial source-space slices through the maximum
activity location by default.

These are source-space activity maps, not anatomical MRI overlays. They are
useful for seeing medial or internal activity that a surface plot can hide.

The MRI slice image ends with:

```text
_mri_slices.png
```

It shows the same source activity overlaid on a brain-only mask from the
standard DIPFIT template MRI. The brain outline comes from the matching
FieldTrip/DIPFIT segmentation file, not from the sparse source grid, so the
edges look more like a natural brain-extracted MRI. The default figure is a compact `2 x 3`
neuroimaging-style montage: sagittal, coronal, and axial slices for the
strongest source location on the first row, and a second spatially distinct
source peak on the second row. The activity is smoothed onto the MRI slice so
the plot does not show jagged source-grid interpolation lines.
This gives a more anatomical view while keeping the toolbox independent of
EEGLAB after the template files have been copied into this folder.

By default, MRI slice rows are chosen from the strongest source activity above
`z = -35 mm`. This avoids ugly inferior template slices that can include
face/neck anatomy. If you want the raw maximum regardless of location, use:

```matlab
'AutoSliceMinZ', -Inf
```

To choose your own slices, add coordinates in millimeters:

```matlab
outputs = doPlotsLoreta(exampleData, ...
    'SliceX', -10, ...
    'SliceY', 20, ...
    'SliceZ', 40);
```

Or set all three at once:

```matlab
outputs = doContrastLoreta(exampleData, ...
    'SliceXYZ', [-10 20 40]);
```

You can also provide more than one MRI slice row:

```matlab
doPlotLoretaMriSlices(outputs.matFile, ...
    'SliceXYZ', [-10 20 40; 20 -30 30], ...
    'OutputPng', 'outputs/custom_mri_slices.png');
```

To turn the automatic slice figure off:

```matlab
'MakeSlicePlot', false
```

To turn the automatic MRI slice figure off:

```matlab
'MakeMriSlicePlot', false
```

By default, figures stay open in MATLAB after they are saved. For batch
processing, you can close figures automatically:

```matlab
'CloseFigures', true
```

You can also plot slices later from a saved output `.mat` file:

```matlab
doPlotLoretaSlices('outputs/rewp_sub01_cond01_source.mat', ...
    'SliceXYZ', [-10 20 40], ...
    'OutputPng', 'outputs/custom_slices.png');
```

For an MRI overlay:

```matlab
doPlotLoretaMriSlices('outputs/rewp_sub01_cond01_source.mat', ...
    'SliceXYZ', [-10 20 40], ...
    'OutputPng', 'outputs/custom_mri_slices.png');
```

To make the MRI overlay smoother or sharper:

```matlab
'OverlaySigmaMm', 12   % smoother
'OverlaySigmaMm', 6    % sharper
```

To change how many automatic MRI peak rows are shown:

```matlab
'NumSliceRows', 1      % one sagittal/coronal/axial row
'NumSliceRows', 2      % default montage
```

### Region CSV

The CSV file ends with:

```text
_regions.csv
```

It ranks brain regions from largest to smallest activity.

Columns include:

- `rank`
- `atlas`
- `label`
- `region`
- `meanActivity`
- `maxActivity`
- `sumActivity`
- `nVertices`

By default, the region ranking uses the Brainstorm `DKT` atlas.

### Group-Level MAT File

The `.mat` file ends with:

```text
_source.mat
```

It contains:

- `sourceERP`
  - the source-localized data
- `sourceValues`
  - the source values at the plotted time point
- `regionRanking`
  - the ranked brain regions
- `metadata`
  - subject, condition, time point, and settings
- `inverse`
  - the inverse model information, including source coordinates

Use this `.mat` file later for group-level analysis.

## Running Different Subjects Or Conditions

Change these numbers:

```matlab
'SubjectIdx', 1
'ConditionIdx', 1
```

## Processing All Subjects For Condition 1

Use:

```matlab
processCondition1AllSubjects
```

At the top of that script, set:

```matlab
timeIndex = 317;
```

That same time point will be used for every subject.

## Recommended Contrast Workflow

For most experiments, the recommended analysis is a condition contrast:

```text
condition 1 - condition 2
```

The code does not care what the conditions mean. They could be wins/losses,
targets/standards, congruent/incongruent, or anything else. Arrange your EEG
matrix so condition 1 is the positive condition and condition 2 is the
negative condition.

Run one subject:

```matlab
exampleContrast
```

Run all subjects:

```matlab
processContrastAllSubjects
```

Run the group-level contrast:

```matlab
exampleGroup
```

Contrast colors:

- warm colors mean `condition 1 > condition 2`
- cool colors mean `condition 2 > condition 1`

## Group-Level Significant-Region Map

After processing all subjects, run:

```matlab
exampleGroup
```

This reads the subject-level `*_contrast.mat` files and creates:

- `*_significant_regions.png`
  - six-view brain map showing only statistically significant regions
- `*_uncorrected_pmap.png`
  - exploratory six-view map colored by uncorrected evidence strength
- `*_group_slices.png`
  - group-average contrast slices through the source grid
- `*_group_mri_slices.png`
  - group-average contrast slices overlaid on the standard template MRI
- `*_all_regions.csv`
  - all tested regions and their statistics, including significance flags
- `*_group.mat`
  - all group-level results, including nonsignificant regions

The current group contrast example uses:

```matlab
inputDir = fullfile('outputs', 'contrast_cond01_minus_cond02_time317');
alpha = 0.05;
correction = 'fdr';
tail = 'both';
```

The group map only colors regions that are statistically significant after
correction. If no regions are significant, the group brain map will be gray.

The uncorrected p-map is descriptive. For contrasts, warm colors mean
`condition 1 > condition 2`, cool colors mean `condition 2 > condition 1`,
and stronger color means smaller uncorrected p value.

The group `.mat` file also stores source-level group statistics. You can
make additional group slice figures from it:

```matlab
doPlotLoretaSlices(outputs.matFile, 'MapType', 't');
doPlotLoretaSlices(outputs.matFile, 'MapType', 'p');
doPlotLoretaSlices(outputs.matFile, 'MapType', 'significant');
```

The same `MapType` options work for MRI overlays:

```matlab
doPlotLoretaMriSlices(outputs.matFile, 'MapType', 't');
doPlotLoretaMriSlices(outputs.matFile, 'MapType', 'p');
doPlotLoretaMriSlices(outputs.matFile, 'MapType', 'significant');
```

For example, subject 4, condition 2 from your own full group file:

```matlab
outputs = doPlotsLoreta('yourGroupData.mat', ...
    'DataVariable', 'yourDataVariable', ...
    'SubjectIdx', 4, ...
    'ConditionIdx', 2, ...
    'OutputDir', 'outputs', ...
    'OutputPrefix', 'sub04_cond02');
```

## Required Files

These files are needed for the high-level function:

- `exampleSubjectData.mat`
  - small one-subject example data file
- `matlocs.mat`
  - channel labels and channel-location information
- `leadfield.mat`
  - the forward model used for source localization
- `templates/cortex/brainstorm_icbm152_cortex_pial_low.mat`
  - the cortex model and brain-region atlases used for plotting and CSV
- `templates/dipfit/standard_BEM/standard_mri.mat`
  - standard template MRI used for MRI slice overlays
- `templates/dipfit/standard_BEM/standard_seg.mat`
  - matching template segmentation used to make brain-only MRI slice overlays

The file `leadfield.mat` has already been generated for the current 63-channel
setup. You do not need to regenerate it unless the channel montage, channel
order, head model, or source grid changes.

## Channel Counts And Leadfields

The EEG data and leadfield must match exactly:

```text
number of EEG channels = number of leadfield rows
```

The channel order must also match. A 63-channel leadfield should not be used
with 31- or 32-channel data.

For 31- or 32-channel data, first create a matching channel-location file,
for example:

```text
matlocs_32chan.mat
```

That file must contain a `chanlocs` variable with the channels in the same
order as the rows of your EEG matrix.

Then build a matching leadfield once:

```matlab
doMakeLeadfield('ChanlocsFile', 'matlocs_32chan.mat', ...
    'OutputFile', 'leadfield_32chan.mat', ...
    'FieldTripPath', '/path/to/fieldtrip');
```

After that, use the matching files when running the analysis:

```matlab
outputs = doPlotsLoreta(data32, ...
    'ChanlocsFile', 'matlocs_32chan.mat', ...
    'LeadfieldFile', 'leadfield_32chan.mat');
```

The same rule applies to contrasts:

```matlab
outputs = doContrastLoreta(data32, ...
    'ChanlocsFile', 'matlocs_32chan.mat', ...
    'LeadfieldFile', 'leadfield_32chan.mat');
```

If the data and leadfield do not match, the toolbox stops with a clear
channel-count or channel-order error.

## If leadfield.mat Is Missing

If `leadfield.mat` is missing, create it once with:

```matlab
doMakeLeadfield
```

This uses local copies of EEGLAB/DIPFIT template files stored in:

```text
templates/dipfit
```

FieldTrip is needed only for this lead-field creation step. Once
`leadfield.mat` exists, the main sLORETA plotting function uses the saved
lead field.

## Advanced Options

Change the atlas:

```matlab
'AtlasName', 'Desikan-Killiany'
```

Other available atlases depend on the Brainstorm cortex template, and include
options such as `DKT`, `Destrieux`, `Desikan-Killiany`, and `Brodmann`.

Change the plot threshold:

```matlab
'ThresholdPercentile', 84
```

Higher values show less activity. Lower values show more activity.

Change MRI slice behavior:

```matlab
'AutoSliceMinZ', -35
'OverlaySigmaMm', 10
'NumSliceRows', 2
'MinPeakDistanceMm', 35
'BrainMaskLabel', 3
'BrainMaskDilateVoxels', 2
```

`AutoSliceMinZ` controls the lowest z-coordinate allowed when the MRI slice
is chosen automatically. `OverlaySigmaMm` controls smoothing of the source
activity on MRI slices. `NumSliceRows` controls how many sagittal/coronal/
axial peak rows are shown. `MinPeakDistanceMm` controls how far apart
automatic peaks must be. `BrainMaskLabel` selects the brain compartment in
`standard_seg.mat`. `BrainMaskDilateVoxels` keeps a small anatomical rim
around the brain mask so the MRI edge does not look artificially chopped.

Restrict the plot to the ACC bounding box:

```matlab
'Region', 'acc'
```

## Lower-Level Functions

Most users do not need these directly, but they are available:

- `sloreta_make_inverse.m`
  - builds the inverse operator from the lead field
- `sloreta_apply.m`
  - applies the inverse operator to EEG/ERP data
- `plot_sloreta_cortex.m`
  - makes the six-view cortex plot
- `doPlotLoretaSlices.m`
  - makes sagittal/coronal/axial internal source-slice plots
- `doPlotLoretaMriSlices.m`
  - overlays sagittal/coronal/axial source maps on a standard template MRI
- `doMakeLeadfield.m`
  - builds a leadfield for a specific channel montage
- `export_sloreta_nifti.m`
  - exports a NIfTI volume for neuroimaging viewers
