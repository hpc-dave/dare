/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef FILEFORMATS_XML_H_
#define FILEFORMATS_XML_H_

#include <iostream>
#include <string>
#include <vector>

namespace dare::ff {

struct _XMLAttribute {
    _XMLAttribute(std::string _name, std::string _value) : name(_name), value(_value) {}
    _XMLAttribute() : name(""), value("") {}
    std::string name;
    std::string value;
};

class XMLNode {
public:
    XMLNode();
    explicit XMLNode(const std::string& _tag);
    XMLNode(const std::string& _tag, const std::string& _data);
    XMLNode(const std::string& _tag, const XMLNode& _nodes);
    XMLNode(const std::string& _tag, const std::vector<XMLNode>& _nodes);
    XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att);
    XMLNode(const std::string& _tag, const _XMLAttribute& _att);
    XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att, const std::vector<XMLNode>& _nodes);
    XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att, const std::string& data);

    XMLNode& operator<<(const XMLNode& node);
    XMLNode& AddChild(const XMLNode& node);

    void AppendToOutput(const std::size_t indent_level, std::ostream& os);

    static void SetIndentLength(std::size_t length);
    static std::size_t GetIndentLength();

private:
    static std::size_t indent_length;
    std::string tag;
    std::vector<_XMLAttribute> attributes;
    std::string data;
    std::vector<XMLNode> children;
};

class XML {
public:
    XML();
    explicit XML(const XMLNode& node);

    XML& AddNode(const XMLNode& node);
    XML& operator<<(const XMLNode& node);
    XML& operator<<(const XML& other);
    XML& AddXMLObject(const XML& other);

    void WriteToFile(const std::string& file_name, std::ios::openmode mode = std::ios::out);

private:
    std::vector<XMLNode> nodes;
};

}  // namespace dare::ff

#endif  // FILEFORMATS_XML_H_
