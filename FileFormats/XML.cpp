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

#include <fstream>
#include "XML.h"

namespace dare::ff{

std::size_t XMLNode::indent_length = 2;

XMLNode::XMLNode(){}

XMLNode::XMLNode(const std::string& _tag)
    : tag(_tag) {
}

XMLNode::XMLNode(const std::string& _tag, const std::string& _data)
    : tag(_tag), data(_data) {
}

XMLNode::XMLNode(const std::string& _tag, const XMLNode& _node) : XMLNode(_tag, std::vector<XMLNode>({_node})) {}

XMLNode::XMLNode(const std::string& _tag, const std::vector<XMLNode>& _nodes)
    : tag(_tag), children(_nodes) {
}

XMLNode::XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att) : tag(_tag), attributes(_att){
}

XMLNode::XMLNode(const std::string& _tag, const _XMLAttribute& _att)
    : XMLNode(_tag, std::vector<_XMLAttribute>({_att})) {
}

XMLNode::XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att, const std::vector<XMLNode>& _nodes)
    : tag(_tag), attributes(_att), children(_nodes) {}

XMLNode::XMLNode(const std::string& _tag, const std::vector<_XMLAttribute>& _att, const std::string& _data)
    : tag(_tag), attributes(_att), data(_data){
}

void XMLNode::SetIndentLength(std::size_t length){
    indent_length = length;
}

std::size_t XMLNode::GetIndentLength(){
    return indent_length;
}

XMLNode& XMLNode::operator<<(const XMLNode& node){
    return AddChild(node);
}

XMLNode& XMLNode::AddChild(const XMLNode& node){
    children.push_back(node);
    return *this;
}

void XMLNode::AppendToOutput(const std::size_t indent_level, std::ostream& os){
    std::string indent(indent_level * indent_length, ' ');
    const std::string newl("\n");

    bool is_empty{data.empty() && children.empty()};

    os << indent << '<' << tag;

    if(is_empty){
        os << " />" << newl;
        return;
    }

    if(!attributes.empty()){
        for(const auto& a : attributes){
            os << ' ' << a.name << '=' << '\"' << a.value << '\"';
        }
    }

    os << ">";

    if(!data.empty()){
        os << data;
    } else {
        os << newl;
        for (auto& child : children) {
            child.AppendToOutput(indent_level + indent_length, os);
        }
    }

    os << "</" << tag << '>' << newl;
}

XML::XML(){
}

XML::XML(const XMLNode& node) : nodes(std::vector<XMLNode>({node})){
}

XML& XML::AddNode(const XMLNode& node){
    nodes.push_back(node);
    return *this;
}

XML& XML::operator<<(const XMLNode& node){
    return AddNode(node);
}

XML& XML::operator<<(const XML& other){
    return AddXMLObject(other);
}

XML& XML::AddXMLObject(const XML& other){
    nodes.insert(nodes.end(), other.nodes.begin(), other.nodes.end());
    return *this;
}

void XML::WriteToFile(const std::string& file_name, std::ios::openmode mode){
    std::ofstream ofs(file_name, std::ios::out);
    if(!ofs){
        std::cerr << "Couldn't open file " << file_name << " aborting writing" << std::endl;
    }

    for(auto& n : nodes){
        n.AppendToOutput(0, ofs);
    }

    ofs.close();
}

} // namespace dare::ff
