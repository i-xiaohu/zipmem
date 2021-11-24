cmake_minimum_required(VERSION 3.10)
project(zipmem C)

set(CMAKE_C_STANDARD 99)
set(exe zipmem)
add_compile_options(-g -O2)

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
        FM_index/bwt.h
        FM_index/bwt_gen.c
        FM_index/is.c
        FM_index/QSufSort.c
        FM_index/QSufSort.h
        FM_index/rle.c
        FM_index/rle.h
        FM_index/rope.c
        FM_index/rope.h)

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
        bwamem/bwamem_extra.c
        bwamem/bwamem_pair.c)

add_executable(bwamem
        ${cstl}
        ${FM_index}
        ${bwalib}
        ${bwamem}
        bwamem/main_mem.c
        read_input.c
        read_input.h)
target_link_libraries(bwamem -lm -lz -lpthread -lrt)

add_executable(${exe}
        ${cstl}
        ${FM_index}
        ${bwalib}
        ${bwamem}
        zipmem.c
        zipmem.h
        read_input.c
        read_input.h
        main.c)
target_link_libraries(${exe} -lm -lz -lpthread -lrt)

add_executable(check
        ${cstl}
        ${FM_index}
        ${bwalib}
        samop.c
        samop.h
        check_seed.c)
target_link_libraries(check -lz -lpthread -lrt)

add_executable(bwaidx
        ${FM_index}
        ${bwalib}
        FM_index/main.c)
target_link_libraries(bwaidx -lz -lrt)