# app.py
# Streamlit UI to generate a unique IFC model from three digits (X,Y,Z in 0..9)

import math, sys, io, tempfile, os
import streamlit as st
import ifcopenshell
from ifcopenshell.api import run

st.set_page_config(page_title="XYZ → IFC", page_icon="🏗️", layout="centered")
st.title("🏗️ DiTWP-UE Hausübung IFC Modell Generator")

st.write("Eingabe der letzten drei Nummern der Matrikelnummer, z.B. 123456XYZ → XYZ, ")

# ── sliders (0..9) ─────────────────────────────────────────
colA, colB, colC = st.columns(3)
with colA: Xd = st.slider("X", 0, 9, 9)
with colB: Yd = st.slider("Y", 0, 9, 9)
with colC: Zd = st.slider("Z", 0, 9, 9)

student_tag = st.text_input("Eingabe der kompletten Matrikelnummer + Enter", "123456789")

# ── mapping digits → mm parameters ─────────────────────────
def digits_to_params(xd:int, yd:int, zd:int):
    # slab size
    X_mm = 4000.0 + xd * 600.0        # 4.0 m … 9.4 m
    Y_mm = 2500.0 + yd * 350.0        # 2.5 m … 5.65 m
    Z_mm = 180.0  + zd * 30.0         # 180 … 450 (also column/beam base)

    # columns (3) along X inside slab: at margin Z, mid, and (X - Z), centered in Y
    col_side_mm   = Z_mm
    col_height_mm = max(2400.0, 10.0 * Z_mm)
    x_positions   = [Z_mm, 0.5 * X_mm, X_mm - Z_mm]
    y_center      = 0.5 * Y_mm

    # beam section equals column size (depth 1.5× just for variety)
    beam_b_mm = Z_mm
    beam_h_mm = Z_mm * 1.5

    return {
        "SLAB_LENGTH_MM": X_mm,
        "SLAB_WIDTH_MM":  Y_mm,
        "SLAB_THICK_MM":  Z_mm,
        "Z_OFFSETS_MM":   [0.0],
        "COLUMNS": [
            {"name":"C1","x":x_positions[0],"y":y_center,"height":col_height_mm,"b":col_side_mm,"h":col_side_mm},
            {"name":"C2","x":x_positions[1],"y":y_center,"height":col_height_mm,"b":col_side_mm,"h":col_side_mm},
            {"name":"C3","x":x_positions[2],"y":y_center,"height":col_height_mm,"b":col_side_mm,"h":col_side_mm},
        ],
        "BEAM_B_MM": beam_b_mm,
        "BEAM_H_MM": beam_h_mm,
    }

P = digits_to_params(Xd, Yd, Zd)

st.subheader("Geometrie (mm)")
st.write(f"- Bodenplatte **{int(P['SLAB_LENGTH_MM'])} × {int(P['SLAB_WIDTH_MM'])} × {int(P['SLAB_THICK_MM'])}** mm")
st.write(f"- Stützen x = {int(P['COLUMNS'][0]['x'])}, {int(P['COLUMNS'][1]['x'])}, {int(P['COLUMNS'][2]['x'])}; "
         f"y = {int(P['COLUMNS'][0]['y'])}; B/H = {int(P['COLUMNS'][0]['b'])}×{int(P['COLUMNS'][0]['h'])} mm; "
         f"Höhe = {int(P['COLUMNS'][0]['height'])} mm")
st.write(f"- Balken **{int(P['BEAM_B_MM'])} × {int(P['BEAM_H_MM'])}** mm")

# ───────────────── IFC builder (geometry only) ───────────────────────
def build_ifc_bytes(P, project_name):
    mm2m = 1.0 / 1000.0
    length_m = P["SLAB_LENGTH_MM"] * mm2m
    width_m  = P["SLAB_WIDTH_MM"]  * mm2m
    thick_m  = P["SLAB_THICK_MM"]  * mm2m

    # column base = slab top
    SUPPORT_TOP_M = P["SLAB_THICK_MM"] * mm2m

    model   = run("project.create_file", version="IFC4")
    project = run("root.create_entity", model, ifc_class="IfcProject", name=project_name)
    run("unit.assign_unit", model)

    site     = run("root.create_entity", model, ifc_class="IfcSite",            name="Site")
    building = run("root.create_entity", model, ifc_class="IfcBuilding",        name="Building")
    storey   = run("root.create_entity", model, ifc_class="IfcBuildingStorey",  name="Ground Floor")
    run("aggregate.assign_object", model, relating_object=project,  products=[site])
    run("aggregate.assign_object", model, relating_object=site,     products=[building])
    run("aggregate.assign_object", model, relating_object=building, products=[storey])

    model3d = run("context.add_context", model, context_type="Model")
    body_ctx = run("context.add_context", model,
                   context_type="Model", context_identifier="Body",
                   target_view="MODEL_VIEW", parent=model3d)

    # helpers
    def make_rect_profile(b_m: float, h_m: float, name: str):
        x = b_m * 0.5; y = h_m * 0.5
        pts = [(-x,-y),(x,-y),(x,y),(-x,y),(-x,-y)]
        return run("profile.add_arbitrary_profile", model, profile=pts, name=name)

    def placement_matrix(origin, x_axis, y_axis, z_axis):
        return (
            (x_axis[0], y_axis[0], z_axis[0], origin[0]),
            (x_axis[1], y_axis[1], z_axis[1], origin[1]),
            (x_axis[2], y_axis[2], z_axis[2], origin[2]),
            (0.0, 0.0, 0.0, 1.0),
        )

    def z_aligned_basis(dir_vec, up_hint=(0.0,0.0,1.0)):
        def norm(v):
            l = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) or 1.0
            return (v[0]/l, v[1]/l, v[2]/l)
        def cross(a,b):
            return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        def dot(a,b):
            return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
        z = norm(dir_vec)
        up = up_hint
        if abs(dot(z, up)) > 0.99:
            up = (0.0, 1.0, 0.0)
        x = norm(cross(up, z))
        y = cross(z, x)
        return x, y, z

    # slabs
    slab_outline = [(0.0,0.0),(length_m,0.0),(length_m,width_m),(0.0,width_m),(0.0,0.0)]
    slab_profile = run("profile.add_arbitrary_profile", model, profile=slab_outline, name="SlabProfile")
    for i, z_mm in enumerate(P["Z_OFFSETS_MM"]):
        rep = run("geometry.add_profile_representation", model, context=body_ctx,
                  profile=slab_profile, depth=thick_m)
        slab = run("root.create_entity", model, ifc_class="IfcSlab", name=f"Slab_{i+1}")
        run("geometry.assign_representation", model, product=slab, representation=rep)
        run("geometry.edit_object_placement", model, product=slab,
            matrix=((1,0,0,0),(0,1,0,0),(0,0,1,z_mm*mm2m),(0,0,0,1)))
        run("spatial.assign_container", model, relating_structure=storey, products=[slab])

    # columns
    def col_top_xyz(col):
        return (col["x"]*mm2m, col["y"]*mm2m, SUPPORT_TOP_M + col["height"]*mm2m)

    for col in P["COLUMNS"]:
        x_m, y_m = col["x"]*mm2m, col["y"]*mm2m
        h_m      = col["height"]*mm2m
        b_m      = col["b"]*mm2m
        d_m      = col["h"]*mm2m

        c_prof = make_rect_profile(b_m, d_m, name=f"{col['name']}_Profile")
        c_rep  = run("geometry.add_profile_representation", model, context=body_ctx,
                     profile=c_prof, depth=h_m)

        column = run("root.create_entity", model, ifc_class="IfcColumn", name=col["name"])
        run("geometry.assign_representation", model, product=column, representation=c_rep)
        run("geometry.edit_object_placement", model, product=column,
            matrix=((1,0,0,x_m),(0,1,0,y_m),(0,0,1,SUPPORT_TOP_M),(0,0,0,1)))
        run("spatial.assign_container", model, relating_structure=storey, products=[column])

    # beams (center-to-center, bottom flush with column tops)
    beam_b_m = P["BEAM_B_MM"] * mm2m
    beam_h_m = P["BEAM_H_MM"] * mm2m
    beam_profile = make_rect_profile(beam_b_m, beam_h_m, name="BeamProfile")
    half_beam_h_m = 0.5 * beam_h_m

    cols = P["COLUMNS"]
    for i in range(len(cols)-1):
        a, b = cols[i], cols[i+1]
        A = col_top_xyz(a)
        B = col_top_xyz(b)

        dir_vec  = (B[0]-A[0], B[1]-A[1], B[2]-A[2])
        span_len = math.sqrt(dir_vec[0]**2 + dir_vec[1]**2 + dir_vec[2]**2)
        if span_len < 1e-9:
            continue

        rep = run("geometry.add_profile_representation", model, context=body_ctx,
                  profile=beam_profile, depth=span_len)
        beam = run("root.create_entity", model, ifc_class="IfcBeam", name=f"Beam_{i+1}")
        run("geometry.assign_representation", model, product=beam, representation=rep)

        x_ax, y_ax, z_ax = z_aligned_basis(dir_vec)
        origin = (A[0] + y_ax[0]*half_beam_h_m,
                  A[1] + y_ax[1]*half_beam_h_m,
                  A[2] + y_ax[2]*half_beam_h_m)
        run("geometry.edit_object_placement", model, product=beam,
            matrix=placement_matrix(origin, x_ax, y_ax, z_ax))
        run("spatial.assign_container", model, relating_structure=storey, products=[beam])

    # return IFC as bytes
    with tempfile.NamedTemporaryFile(delete=False, suffix=".ifc") as tmp:
        model.write(tmp.name)
        tmp.flush()
        with open(tmp.name, "rb") as f: data = f.read()
    os.unlink(tmp.name)
    return data

# ── generate & download ───────────────────────────────────
outfile = f"ifc_{Xd}{Yd}{Zd}.ifc"
proj_name = f"DiTWP IFC [{Xd}{Yd}{Zd}]" + (f" - {student_tag}" if student_tag else "")

if st.button("Erstellen IFC"):
    data = build_ifc_bytes(P, proj_name)
    st.success("IFC erstellt.")
    st.download_button("⬇️ IFC herunterladen", data=data, file_name=outfile, mime="application/octet-stream")
