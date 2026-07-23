// GlLoader.cpp -- see GlLoader.h.

#include "app/gui/GlLoader.h"

#include <GLFW/glfw3.h>

#include <type_traits>

namespace glx {

GLuint (*CreateShader)(GLenum) = nullptr;
void (*ShaderSource)(GLuint, GLsizei, const char* const*, const GLint*) = nullptr;
void (*CompileShader)(GLuint) = nullptr;
void (*GetShaderiv)(GLuint, GLenum, GLint*) = nullptr;
void (*GetShaderInfoLog)(GLuint, GLsizei, GLsizei*, char*) = nullptr;
GLuint (*CreateProgram)() = nullptr;
void (*AttachShader)(GLuint, GLuint) = nullptr;
void (*LinkProgram)(GLuint) = nullptr;
void (*GetProgramiv)(GLuint, GLenum, GLint*) = nullptr;
void (*GetProgramInfoLog)(GLuint, GLsizei, GLsizei*, char*) = nullptr;
void (*DeleteShader)(GLuint) = nullptr;
void (*DeleteProgram)(GLuint) = nullptr;
void (*UseProgram)(GLuint) = nullptr;
GLint (*GetUniformLocation)(GLuint, const char*) = nullptr;
void (*UniformMatrix4fv)(GLint, GLsizei, GLboolean, const GLfloat*) = nullptr;
void (*Uniform1f)(GLint, GLfloat) = nullptr;
void (*Uniform1i)(GLint, GLint) = nullptr;
void (*Uniform2f)(GLint, GLfloat, GLfloat) = nullptr;
void (*Uniform4f)(GLint, GLfloat, GLfloat, GLfloat, GLfloat) = nullptr;
void (*GenBuffers)(GLsizei, GLuint*) = nullptr;
void (*BindBuffer)(GLenum, GLuint) = nullptr;
void (*BufferData)(GLenum, glSizeiptr, const void*, GLenum) = nullptr;
void (*DeleteBuffers)(GLsizei, const GLuint*) = nullptr;
void (*GenVertexArrays)(GLsizei, GLuint*) = nullptr;
void (*BindVertexArray)(GLuint) = nullptr;
void (*DeleteVertexArrays)(GLsizei, const GLuint*) = nullptr;
void (*BindAttribLocation)(GLuint, GLuint, const char*) = nullptr;
void (*EnableVertexAttribArray)(GLuint) = nullptr;
void (*VertexAttribPointer)(GLuint, GLint, GLenum, GLboolean, GLsizei, const void*) = nullptr;
void (*GenFramebuffers)(GLsizei, GLuint*) = nullptr;
void (*BindFramebuffer)(GLenum, GLuint) = nullptr;
void (*FramebufferTexture2D)(GLenum, GLenum, GLenum, GLuint, GLint) = nullptr;
void (*FramebufferRenderbuffer)(GLenum, GLenum, GLenum, GLuint) = nullptr;
GLenum (*CheckFramebufferStatus)(GLenum) = nullptr;
void (*DeleteFramebuffers)(GLsizei, const GLuint*) = nullptr;
void (*GenRenderbuffers)(GLsizei, GLuint*) = nullptr;
void (*BindRenderbuffer)(GLenum, GLuint) = nullptr;
void (*RenderbufferStorage)(GLenum, GLenum, GLsizei, GLsizei) = nullptr;
void (*DeleteRenderbuffers)(GLsizei, const GLuint*) = nullptr;

bool init() {
    static bool done = false, ok = false;
    if (done) return ok;
    done = true;
    bool all = true;
    auto load = [&](auto& fn, const char* name) {
        fn = reinterpret_cast<std::remove_reference_t<decltype(fn)>>(
            glfwGetProcAddress(name));
        if (!fn) all = false;
    };
    load(CreateShader, "glCreateShader");
    load(ShaderSource, "glShaderSource");
    load(CompileShader, "glCompileShader");
    load(GetShaderiv, "glGetShaderiv");
    load(GetShaderInfoLog, "glGetShaderInfoLog");
    load(CreateProgram, "glCreateProgram");
    load(AttachShader, "glAttachShader");
    load(LinkProgram, "glLinkProgram");
    load(GetProgramiv, "glGetProgramiv");
    load(GetProgramInfoLog, "glGetProgramInfoLog");
    load(DeleteShader, "glDeleteShader");
    load(DeleteProgram, "glDeleteProgram");
    load(UseProgram, "glUseProgram");
    load(GetUniformLocation, "glGetUniformLocation");
    load(UniformMatrix4fv, "glUniformMatrix4fv");
    load(Uniform1f, "glUniform1f");
    load(Uniform1i, "glUniform1i");
    load(Uniform2f, "glUniform2f");
    load(Uniform4f, "glUniform4f");
    load(GenBuffers, "glGenBuffers");
    load(BindBuffer, "glBindBuffer");
    load(BufferData, "glBufferData");
    load(DeleteBuffers, "glDeleteBuffers");
    load(GenVertexArrays, "glGenVertexArrays");
    load(BindVertexArray, "glBindVertexArray");
    load(DeleteVertexArrays, "glDeleteVertexArrays");
    load(BindAttribLocation, "glBindAttribLocation");
    load(EnableVertexAttribArray, "glEnableVertexAttribArray");
    load(VertexAttribPointer, "glVertexAttribPointer");
    load(GenFramebuffers, "glGenFramebuffers");
    load(BindFramebuffer, "glBindFramebuffer");
    load(FramebufferTexture2D, "glFramebufferTexture2D");
    load(FramebufferRenderbuffer, "glFramebufferRenderbuffer");
    load(CheckFramebufferStatus, "glCheckFramebufferStatus");
    load(DeleteFramebuffers, "glDeleteFramebuffers");
    load(GenRenderbuffers, "glGenRenderbuffers");
    load(BindRenderbuffer, "glBindRenderbuffer");
    load(RenderbufferStorage, "glRenderbufferStorage");
    load(DeleteRenderbuffers, "glDeleteRenderbuffers");
    ok = all;
    return ok;
}

}  // namespace glx
