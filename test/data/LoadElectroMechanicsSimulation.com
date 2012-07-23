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

# Set background colour to white.
gfx modify window 1 background colour 1 1 1 texture none;

# Make a colour bar
gfx create colour_bar spectrum default number_format %2.f label_material black extend_length 0

# Add the colour bar to the scene
gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre -2.6,0,0 select_on material default selected_material default normalised_window_fit_left;


