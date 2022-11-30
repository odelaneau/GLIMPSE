#ifndef __MAKE_UNIQUE_HPP__
#define __MAKE_UNIQUE_HPP__

#include <memory>

#if __cplusplus < 201402L

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

#else

using std::make_unique;

#endif

#endif /* __MAKE_UNIQUE_HPP__ */