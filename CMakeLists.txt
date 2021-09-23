cmake_minimum_required(VERSION 3.10)
project(zip_seeding C)

set(CMAKE_C_STANDARD 99)
set(exe zipmem)
add_compile_options(-O2)

set(cstl
        cstl/kbtree.h
        cstl/khash.h
        cstl/ksort.h
        cstl/kstring.c
        cstl/kstring.h
        cstl/kthread.c
        cstl/kthread.h
        cstl/kvec.h)

set(FM_index
        FM_index/bntseq.c
        FM_index/bntseq.h
        FM_index/bwt.c
        FM_index/bwt.h)

set(bwalib
        bwalib/bwa.c
        bwalib/bwa.h
        bwalib/bwashm.c
        bwalib/kopen.c
        bwalib/kopen.h
        bwalib/kseq.h
        bwalib/ksw.c
        bwalib/ksw.h
        bwalib/utils.c
        bwalib/utils.h)

set(bwamem
        bwamem/bwamem.c
        bwamem/bwamem.h
        )

add_executable(bwamem
        ${cstl}
        ${FM_index}
        ${bwalib}
        ${bwamem}
        bwamem/main_mem.c)
target_link_libraries(bwamem -lm -lz -lpthread -lrt)

add_executable(${exe}
        ${cstl}
        ${FM_index}
        ${bwalib}
        ${bwamem}
        zipmem.c
        zipmem.h
        main.c)
target_link_libraries(${exe} -lm -lz -lpthread -lrt)