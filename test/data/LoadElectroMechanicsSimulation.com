gfx read node undeformed time -1
for ($i=0; $i<=36; $i++) { 
  gfx read node solution_$i time $i
  gfx read node ../../electrics/cmgui_output/voltage_mechanics_mesh_$i time $i
}
gfx read ele undeformed
gfx define faces egroup solution




gfx modify g_element solution general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element solution surfaces select_on material default data V spectrum default selected_material default_selected render_shaded;

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -90 60 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx cr win
