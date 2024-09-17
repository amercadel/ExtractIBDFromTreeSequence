project('extractIBD', ['c', 'cpp'],
    default_options: ['c_std=c99', 'cpp_std=c++17'])


kastore_inc = include_directories('kastore/c')
kastore = static_library('kastore', 'kastore/c/kastore.c')
kastore_dep = declare_dependency(link_with: kastore, include_directories: kastore_inc)


cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)
tskit_inc = [kastore_inc, include_directories('tskit/c')]
lib_deps = [m_dep, kastore_dep, declare_dependency(include_directories: tskit_inc)]
lib_sources = [
    'tskit/c/tskit/core.c', 'tskit/c/tskit/tables.c', 'tskit/c/tskit/trees.c',
    'tskit/c/tskit/genotypes.c', 'tskit/c/tskit/stats.c', 'tskit/c/tskit/convert.c', 'tskit/c/tskit/haplotype_matching.c'
]

omp = dependency('openmp')



tskit_lib = static_library('tskit',
    sources: lib_sources, dependencies: lib_deps)
tskit_dep = declare_dependency(include_directories:tskit_inc, link_with: tskit_lib)
executable('extractIBD', 'src/main.cpp', 'src/read_rate_map.cpp',
            dependencies: [tskit_dep, omp])

