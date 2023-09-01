using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace physics_equations
{
    public class QaSystem
    {
        string[] context;
        public QaSystem(string[] contextLines)
        {
            context = contextLines;
        }

        class WordCounter
        {

            public Dictionary<string, int> counts = new Dictionary<string, int>();

            public void Increment(string word)
            {
                if (counts.ContainsKey(word))
                {
                    counts[word] += 1;
                }
                else
                {
                    counts[word] = 1;
                }
            }

            public int this[string word]
            {
                get
                {
                    if (counts.ContainsKey(word))
                    {
                        return counts[word];
                    }
                    else
                    {
                        return 0;
                    }
                }
            }
        }

        float JaccardSimilarity(string[] query, string[] document)
        {
            HashSet<string> querySet = new HashSet<string>(query);
            HashSet<string> documentSet = new HashSet<string>(document);

            querySet.IntersectWith(documentSet);

            float intersection = querySet.Count;
            float union = querySet.Count + documentSet.Count - intersection;

            return intersection / union;
        }
        float GetCosineSimilarity(string[] q, string[] c)
        {

            WordCounter question = new WordCounter();
            WordCounter context = new WordCounter();

            // Increment counts for question        
            foreach (string w in q)
            {
                question.Increment(w);
            }

            // Increment counts for context        
            foreach (string w in c)
            {
                context.Increment(w);
            }

            float dotProduct = 0;
            float questionMagnitude = 0;
            float contextMagnitude = 0;

            // Loop through all words seen so far        
            foreach (KeyValuePair<string, int> pair in question.counts)
            {

                string w = pair.Key;
                int qCount = question[w];
                int cCount = context[w];

                dotProduct += qCount * cCount;
                questionMagnitude += qCount * qCount;
                contextMagnitude += cCount * cCount;
            }
            float denominator = questionMagnitude * contextMagnitude;

            // Calculate cosine similarity         
            float cosineSimilarity = dotProduct / denominator;

            return cosineSimilarity;
        }
        float GetSemanticSimilarity(string[] q, string[] c)
        {
            double similarity = 0;

            // Exact word matches 
            foreach (string qw in q)
            {
                if (c.Contains(qw))
                {
                    similarity += 0.8;
                }
            }

            // Number of shared words         
            HashSet<string> sharedWords = new HashSet<string>(q);
            sharedWords.IntersectWith(c);
            similarity += sharedWords.Count * 0.2;

            return (float)similarity;
        }

        string[] deleteWords = new string[]
        {
            "a",
            "an",
            "the",
            "of",
            "and",
            "or",
            "as",
            "by",
            "for",
            "it",
            "in",
            "into",
            "on",
            "to",
            "with",
            "who"
        };
        string[] TokenizeQuestion(string question)
        {

            // Replace punctuation with space  
            question = question.ToLower();
            question = Regex.Replace(question, @"\p{P}", " ");

            // Split on spaces      
            string[] words = question.Split(' ');

            // Delete words in deleteWords array
            List<string> tokens = new List<string>();
            foreach (string word in words)
            {
                if (!deleteWords.Contains(word.ToLower()))
                {
                    tokens.Add(word);
                }
            }

            return tokens.ToArray();

        }
        public string GetAnswerCosine(string question)
        {
            string bestMatch = "";
            float maxSimilarity = 0;

            // Tokenize question       
            string[] questionTokens = TokenizeQuestion(question);

            foreach (string line in context)
            {
                // Tokenize context line
                string[] contextTokens = line.Split(' ');

                // Calculate semantic similarity             
                float similarity = GetCosineSimilarity(questionTokens, contextTokens);

                if (similarity > maxSimilarity)
                {
                    maxSimilarity = similarity;
                    bestMatch = line;
                }
            }

            // Fallback if no good match    
            if (bestMatch == "") bestMatch = "I don't know...";

            return bestMatch;
        }
        public string GetAnswerSemantic(string question)
        {
            string bestMatch = "";
            float maxSimilarity = 0;

            // Tokenize question       
            string[] questionTokens = TokenizeQuestion(question);

            foreach (string line in context)
            {
                // Tokenize context line
                string[] contextTokens = line.Split(' ');

                // Calculate semantic similarity             
                float similarity = GetSemanticSimilarity(questionTokens, contextTokens);

                if (similarity > maxSimilarity)
                {
                    maxSimilarity = similarity;
                    bestMatch = line;
                }
            }

            // Fallback if no good match    
            if (bestMatch == "") bestMatch = "I don't know...";

            return bestMatch;
        }
        public string GetAnswerJaccard(string question)
        {
            string bestMatch = "";
            float maxSimilarity = 0;

            // Tokenize question       
            string[] questionTokens = TokenizeQuestion(question);

            foreach (string line in context)
            {
                // Tokenize context line
                string[] contextTokens = line.Split(' ');

                // Calculate semantic similarity             
                float similarity = JaccardSimilarity(questionTokens, contextTokens);

                if (similarity > maxSimilarity)
                {
                    maxSimilarity = similarity;
                    bestMatch = line;
                }
            }

            // Fallback if no good match    
            if (bestMatch == "") bestMatch = "I don't know...";

            return bestMatch;
        }
        float GetLevenSimilarity(string[] tokens1, string[] tokens2)
        {
            int n = tokens1.Length;
            int m = tokens2.Length;
            int[,] distance = new int[n + 1, m + 1];

            // Initialize distance matrix
            for (int i = 0; i <= n; i++)
            {
                distance[i, 0] = i;
            }
            for (int j = 0; j <= m; j++)
            {
                distance[0, j] = j;
            }

            // Calculate Levenshtein distance
            for (int j = 1; j <= m; j++)
            {
                for (int i = 1; i <= n; i++)
                {
                    int cost = (tokens1[i - 1] == tokens2[j - 1]) ? 0 : 1;
                    distance[i, j] = Math.Min(Math.Min(
                        distance[i - 1, j] + 1,
                        distance[i, j - 1] + 1),
                        distance[i - 1, j - 1] + cost);
                }
            }

            // Calculate similarity score
            int maxLen = Math.Max(n, m);
            float similarity = 1 - (float)distance[n, m] / maxLen;
            return similarity;
        }
        public string GetAnswerLeven(string question)
        {
            string bestMatch = "";
            float maxSimilarity = 0;

            // Tokenize question       
            string[] questionTokens = TokenizeQuestion(question);

            foreach (string line in context)
            {
                // Tokenize context line
                string[] contextTokens = line.Split(' ');

                // Calculate semantic similarity             
                float similarity = GetLevenSimilarity(questionTokens, contextTokens);

                if (similarity > maxSimilarity)
                {
                    maxSimilarity = similarity;
                    bestMatch = line;
                }
            }

            // Fallback if no good match    
            if (bestMatch == "") bestMatch = "I don't know...";

            return bestMatch;
        }

    }
}
