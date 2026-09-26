#pragma once
#include <cstddef>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

struct Guarded {
  unsigned char* base;
  unsigned char* data;
  size_t allocated;
  Guarded(size_t bytes, bool at_end) {
#ifdef _WIN32
    SYSTEM_INFO si{};
    GetSystemInfo(&si);
    const size_t page = si.dwPageSize;
#else
    const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
#endif
    const size_t usable = (bytes + page - 1) / page * page;
    allocated = usable + 2 * page;
#ifdef _WIN32
    base = static_cast<unsigned char*>(VirtualAlloc(nullptr, allocated, MEM_RESERVE, PAGE_NOACCESS));
    if (!base || !VirtualAlloc(base + page, usable, MEM_COMMIT, PAGE_READWRITE))
      throw std::runtime_error("VirtualAlloc");
#else
    base = static_cast<unsigned char*>(mmap(nullptr, allocated, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    if (base == MAP_FAILED || mprotect(base + page, usable, PROT_READ | PROT_WRITE))
      throw std::runtime_error("mmap");
#endif
    data = base + page + (at_end ? usable - bytes : 0);
  }
  ~Guarded() {
#ifdef _WIN32
    VirtualFree(base, 0, MEM_RELEASE);
#else
    munmap(base, allocated);
#endif
  }
};
