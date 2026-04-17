#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

class ConfigFile
{
private:
    std::map<std::string, std::string> contents;

    static std::string trim(const std::string& s)
    {
        std::size_t start = 0;
        while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
            ++start;

        std::size_t end = s.size();
        while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
            --end;

        return s.substr(start, end - start);
    }

    void parseLine(const std::string& line)
    {
        std::string cleaned = line;

        // Supprimer commentaires
        std::size_t commentPos = cleaned.find('#');
        if (commentPos != std::string::npos)
            cleaned = cleaned.substr(0, commentPos);

        cleaned = trim(cleaned);
        if (cleaned.empty())
            return;

        std::size_t eqPos = cleaned.find('=');
        if (eqPos == std::string::npos)
            return;

        std::string key = trim(cleaned.substr(0, eqPos));
        std::string value = trim(cleaned.substr(eqPos + 1));

        if (!key.empty())
            contents[key] = value;
    }

public:
    ConfigFile() = default;

    explicit ConfigFile(const std::string& filename)
    {
        std::ifstream file(filename.c_str());
        if (!file)
            throw std::runtime_error("Impossible d'ouvrir le fichier de configuration : " + filename);

        std::string line;
        while (std::getline(file, line))
        {
            parseLine(line);
        }
    }

    void process(const std::string& arg)
    {
        parseLine(arg);
    }

    template <typename T>
    T get(const std::string& key) const
    {
        auto it = contents.find(key);
        if (it == contents.end())
            throw std::runtime_error("Parametre manquant : " + key);

        std::istringstream iss(it->second);
        T value;
        iss >> value;

        if (iss.fail())
            throw std::runtime_error("Conversion impossible pour la cle : " + key + " avec valeur : " + it->second);

        return value;
    }

    template <typename T>
    T get(const std::string& key, const T& defaultValue) const
    {
        auto it = contents.find(key);
        if (it == contents.end())
            return defaultValue;

        std::istringstream iss(it->second);
        T value;
        iss >> value;

        if (iss.fail())
            throw std::runtime_error("Conversion impossible pour la cle : " + key + " avec valeur : " + it->second);

        return value;
    }

    std::string get(const std::string& key, const std::string& defaultValue) const
    {
        auto it = contents.find(key);
        if (it == contents.end())
            return defaultValue;
        return it->second;
    }

    std::string get(const std::string& key) const
    {
        auto it = contents.find(key);
        if (it == contents.end())
            throw std::runtime_error("Parametre manquant : " + key);
        return it->second;
    }
};

#endif