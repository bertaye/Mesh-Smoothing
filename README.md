# CatmullClark-Sqrt3-Subdivision

## Catmull Clark Algorithm
#### To use it, load mesh with correct method ( use mesh->loadOffTriMesh("mesh.off") for triangulated meshes and use mesh->loadOffQuadMesh("mesh.off") for quad meshes)
#### If you loaded a triangulated mesh, call the mesh->convertQuadbySplit() method to continue.
#### Then use mesh->CatmullClarkSubdiv() for Catmull Clark algorithm. Each call corresponds one iteration.
##### Note: convertQuadbySplit() uses Catmull-Clark split method.
#### Here are some example results:

### Cube:
<p float="center">
  <img src= "https://user-images.githubusercontent.com/39909689/166225194-7e90d75c-a489-44f8-a9a3-fca08c7c5ee8.png" width="150" />
  <img src="https://user-images.githubusercontent.com/39909689/166225198-12bd871f-7cc1-4bac-b6b8-4a86fb9f3139.png" width="150" /> 
  <img src="https://user-images.githubusercontent.com/39909689/166225191-aef315f7-2458-4ce7-9dea-58a5103dd139.png"width="150" />
</p>

### Cube w/triangles converted to quads by split

<p float="center">
  <img src= "https://user-images.githubusercontent.com/39909689/166225880-a7c87f72-e6f6-46d3-b281-a6f8c914f2f9.png" width="150" />
  <img src="https://user-images.githubusercontent.com/39909689/166225971-30838738-4e57-48a9-8ede-f7dcff6eff6d.png" width="150" /> 
  <img src="https://user-images.githubusercontent.com/39909689/166225982-0125870e-bdf9-4edb-8f97-6bffe35a64de.png"width="100" />
</p>
### Man
<p float="center">
  <img src= "https://user-images.githubusercontent.com/39909689/166226374-9e7e7b60-9edf-4350-bb36-f6a357ae0b33.png" width="150" />
  <img src= "https://user-images.githubusercontent.com/39909689/166226855-8723bc4a-3192-4b78-98b2-da39f1aa9207.png" width="150" />

  <img src="https://user-images.githubusercontent.com/39909689/166226377-80f9d7b9-4518-4095-9297-81500bbe62b8.png" width="150" /> 
  <img src="https://user-images.githubusercontent.com/39909689/166226656-c16b7879-4bdd-4fc8-8961-29ca81baf1d0.png"width="150" />
</p>


## Sqrt3 Subdivision Algorithm

#### To use it, load mesh with correct method ( use mesh->loadOffTriMesh("mesh.off") because this algorithm only works for triangulated meshes and I dont provide any conversion from quads to triangles right now)
#### Then use mesh->Sqrt3SubDiv for Sqrt3 Subdivision algorithm. Each call corresponds one iteration.
#### Here are some example results:
### Coffee Cup

<p float="center">
  <img src= "https://user-images.githubusercontent.com/39909689/166227447-c9ec774d-dd88-45be-9768-0051e5e29f2b.png" width="150" />
  <img src= "https://user-images.githubusercontent.com/39909689/166227441-6b5363c5-e6aa-4b6d-9c28-7cc8becae36a.png" width="150" />

  <img src="https://user-images.githubusercontent.com/39909689/166227440-44f19406-2dfe-4989-ad48-37b5425703d9.png" width="150" /> 
  <img src="https://user-images.githubusercontent.com/39909689/166227435-6983ff14-90af-4f65-a275-b37edd6c0100.png"width="150" />
</p>

### Horse

<p float="center">
  <img src= "https://user-images.githubusercontent.com/39909689/166227885-246ec1e7-10f6-491b-b66e-604fa2d0c976.png" width="150" />
  <img src= "https://user-images.githubusercontent.com/39909689/166227893-6def6682-4bfb-4c24-8548-7bb1569fbe27.png" width="150" />

  <img src="https://user-images.githubusercontent.com/39909689/166227892-d815a5da-6223-4d94-9738-493d8a3946f3.png" width="150" /> 
  <img src="https://user-images.githubusercontent.com/39909689/166227890-9822449d-1e46-4bd5-b84f-3646b767ca33.png"width="150" />
</p>


