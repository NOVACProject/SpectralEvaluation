#pragma once

#include <exception>
#include <string>

namespace novac
{

// General exception signalling that something was not found.
class NotFoundException : public std::exception
{
public:
    NotFoundException(std::string msg) : std::exception(), message(msg) {}

    const std::string message;

    virtual char const* what() const noexcept override
    {
        return message.c_str();
    }
};

// General exception signalling that we failed to read/write a file
class FileIoException : public std::exception
{
public:
    FileIoException(std::string msg) : std::exception(), message(msg) {}

    const std::string message;

    virtual char const* what() const noexcept override
    {
        return message.c_str();
    }
};

}
