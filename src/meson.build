sources = files(
  'main.cpp',
  'instance.cpp',
  'solution.cpp',
  'MLP.cpp',
  'construction.cpp',
  'RVND.cpp',
  'GILS.cpp',
  'neighborhood.cpp',
  'perturbation.cpp',
)

incdir = include_directories('../include')

executable(
  meson.project_name(), 
  sources,
  include_directories: [incdir, incdir_deps],
  dependencies: [dependencies, external_dependencies],
  install: true
)

