/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011, 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::wordRe::meta(char c)
{
    return regExp::meta(c);
}


inline bool Foam::wordRe::valid(char c)
{
    return keyType::valid(c);
}


inline bool Foam::wordRe::isPattern(const std::string& str)
{
    return string::meta<regExp>(str);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordRe::wordRe()
:
    word(),
    re_()
{}


inline Foam::wordRe::wordRe(const wordRe& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(wordRe&& str)
:
    word(std::move(static_cast<word&>(str))),
    re_(std::move(str.re_))
{}


inline Foam::wordRe::wordRe(const keyType& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(const char* str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const std::string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const word& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const keyType& str, const compOption opt)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile(opt);
    }
}


inline Foam::wordRe::wordRe(const char* str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const std::string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const word& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRe::isLiteral() const
{
    return !re_.exists();
}


inline bool Foam::wordRe::isPattern() const
{
    return re_.exists();
}


inline bool Foam::wordRe::compile(const compOption opt)
{
    if (opt)
    {
        bool comp = false;

        if (opt & wordRe::REGEX)
        {
            comp = true;
        }
        else if (opt & wordRe::DETECT)
        {
            comp = string::meta<regExp>(*this) || !string::valid<word>(*this);
        }
        else if (opt & wordRe::ICASE)
        {
            comp = true;
        }

        if (comp)
        {
            return re_.set(*this, (opt & wordRe::ICASE));
        }
    }

    // Fall-through behaviour - not a regex
    re_.clear();
    return false;
}


inline bool Foam::wordRe::compile()
{
    return re_.set(*this);
}


inline void Foam::wordRe::uncompile(bool doStrip)
{
    // Only strip when debug is active (potentially costly operation)
    if (re_.clear() && doStrip && word::debug)
    {
        string::stripInvalid<word>(*this);
    }
}


inline void Foam::wordRe::clear()
{
    word::clear();
    re_.clear();
}


inline bool Foam::wordRe::match(const std::string& text, bool literal) const
{
    if (literal || !re_.exists())
    {
        return !compare(text);  // Compare as literal string
    }

    return re_.match(text);     // Match as regex
}


inline Foam::string Foam::wordRe::quotemeta() const
{
    return string::quotemeta<regExp>(*this);
}


inline void Foam::wordRe::set(const std::string& str, const compOption opt)
{
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::set(const char* str, const compOption opt)
{
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::swap(wordRe& str)
{
    // Self-swapping is a no-op
    if (this == &str)
    {
        return;
    }

    word::swap(static_cast<word&>(str));
    re_.swap(str.re_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRe::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::wordRe::operator=(const wordRe& str)
{
    // Self-assignment is a no-op
    if (this == &str)
    {
        return;
    }

    assign(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const word& str)
{
    assign(str);
    re_.clear();
}


inline void Foam::wordRe::operator=(const keyType& str)
{
    assign(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const std::string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const char* str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(wordRe&& str)
{
    // Self-assignment is a no-op
    if (this == &str)
    {
        return;
    }

    clear();
    swap(str);
}


// ************************************************************************* //
