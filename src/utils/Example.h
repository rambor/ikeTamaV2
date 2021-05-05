//
// Created by xos81802 on 17/01/2020.
//

#ifndef IKETAMA_EXAMPLE_H
#define IKETAMA_EXAMPLE_H

#include <string>

class Example {

    const std::string type;

    void printAnchor();
    void printHelical();
    void printManual();

public:
    Example(std::string what);
};


#endif //IKETAMA_EXAMPLE_H
