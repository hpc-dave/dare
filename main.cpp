/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

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

#include "ScopeGuard/ScopeGuard.h"

int main(int argc, char* argv[]) {
    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        if (scope_guard.AmIRoot()) {
            std::cout << "Hello, nice that you are here!" << std::endl
                      << std::endl;
            std::cout << "  _,-'`''-~`) \n"
                      << "(`~_,=========\\ \n"
                      << " |---,___.-.__,\\ \n"
                      << " |        o     \\ ___  _,,,,_     _.--.\n"
                      << "  \\      `^`    /`_.-'~      `~-;`     \\ \n"
                      << "   \\_      _  .'                 `,     |\n"
                      << "    |`-                           \\'__/\n"
                      << "   /                      ,_       \\  `'-.\n"
                      << "  /    .-''~~--.            `'-,   ;_    /\n"
                      << " |              \\               \\  | `''`\n"
                      << " \\__.--'`'-.   /_               |'\n"
                      << "             `'`  `~~~---..,     |\n"
                      << " jgs                         \\ _.-'`-.\n"
                      << "                             \\       \\ \n"
                      << "                              '.     /\n"
                      << "                                `'~'`\n";
            // copyright by Joan G. Stark, taken from https://www.asciiart.eu/animals/bears
        }
    }
    return 0;
}
