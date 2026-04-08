All-sky 1 degree maps of the synchrotron spectral index and curvature can be downloaded here:

https://drive.google.com/file/d/11N8u-Vrd3stFNN_v6UirGjLbcYrJm5GC/view?usp=sharing

Please cite the following paper: 

https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/stag517/8524011

To produce maps: 
1) Run preprocessing notebooks: preprcessing_XXX_XXX.ipynb for North Coarse, North Fine, South Coarse and South Fine
2) Run LSQ notebooks: coarse_parameters.ipynb, coarse_parameters_southern.ipynb, fine_parameters.ipynb, fine_parameters_southern.ipynb
3) Combine 4 maps into a single map for each parameters (beta and Cs): matching_spectral_maps_v2.ipynb

To plot paper results you can use: show_4fit_maps.ipynb and empirical_comp_gsm.ipynb
