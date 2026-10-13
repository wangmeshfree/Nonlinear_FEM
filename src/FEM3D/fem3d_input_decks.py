"""Student-editable settings; component flags are (x,y,z): 0,1,99."""
try:
    from .fem3d_exact_solutions import UniformTension, ManufacturedSine, ManufacturedZeroTraction
except ImportError:
    from fem3d_exact_solutions import UniformTension, ManufacturedSine, ManufacturedZeroTraction

problem_settings={

    'Uniaxial_tension':{
        'model_characteristics':[1.0],
        'material_properties':[1000.,.3], 
        'integration_order':2,
        'mesh_file':'cube/3d_cube_16.k',
        'mesh_files':['cube/3d_cube_2.k','cube/3d_cube_4.k','cube/3d_cube_8.k'],
        'dirichlet':{'XMIN':(1,1,1)}, 
        'dirichlet_values':{'XMIN':(0.,0.,0.)},
        'traction':{'XMAX':(1,0,0)},
        'traction_values':{'XMAX':(100.,0.,0.)},
        'convergence':False,
        'title':'3D uniform tension (affine patch test)',
        'body_force':(0.,0.,0.)},
         

    'Manufactured_cube':{
        'solution':ManufacturedZeroTraction,'model_characteristics':[.01],
        'material_properties':[1000.,.3],
        'integration_order':2,
        'mesh_file':'cube/3d_cube_4.k',
        'mesh_files':['cube/3d_cube_2.k','cube/3d_cube_4.k','cube/3d_cube_8.k','cube/3d_cube_16.k'],
        'dirichlet':{name:(99,99,99) for name in ['XMIN','XMAX']},
        'traction':{},
        'convergence':True,
        'title':'3D manufactured displacement'},
     
    'Cantilever_solid':{
        'material_properties':[1000.,.3],
        'integration_order':2,
        'mesh_file':'cantilever.k',
        'dirichlet':{'XMIN':(1,1,1)},
        'dirichlet_values':{'XMIN':(0.,0.,0.)},
        'traction':{'XMAX':(0,1,0)},
        'traction_values':{'XMAX':(0.,.01,0.)},
        'body_force':(0.,0.,0.),
        'convergence':False,'title':'3D solid cantilever'},
}
