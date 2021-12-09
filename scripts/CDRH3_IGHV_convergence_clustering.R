### CDRH3 total clone count change nTD vs TD (BY IGHV) ----

library("edgeR")# [1] ‘3.26.1’
library("plyr")
library("dplyr") # [1] ‘1.0.4’
library("data.table") # [1] ‘1.13.6’
library("rlist") # [1] ‘0.4.6.1’
library("textshape") # [1] ‘1.7.1’
library("sjmisc") # [1] ‘2.8.7’
library("reshape")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("Biostrings")
library("gdata")
library("stringr")
library("gridExtra")
library("tidyverse")
library("stringdist")
library("plotly")
library("seqinr")

# Reload metadata and VDJtools output file
# Set background settings
IGHV <- "IGHV3-7" # Or other IGHVs of interest
hamming_distance <- 1 # Set Hamming distance

IgH_CDR3_m <- # File containing BCR information (most importantly CDR H3 sequences)
# Remove values with * or _ etc.
IgH_CDR3_m <- IgH_CDR3_m[!grepl("_", IgH_CDR3_m$cdr3aa),]
IgH_CDR3_m <- IgH_CDR3_m[!grepl("\\*", IgH_CDR3_m$cdr3aa),]
# Filter CDRH3 with IGHV
IgH_CDR3_m_filtered <- IgH_CDR3_m[IgH_CDR3_m$v == IGHV,]
CDR3_length_filtered <- IgH_CDR3_m_filtered
matrix_sequences <- CDR3_length_filtered$cdr3aa
# Hamming distance calculation
CDRH_distance_temp <- stringdistmatrix(matrix_sequences, matrix_sequences, method = "hamming")
# Set limit to hamming distance
CDRH_distance <- as.data.frame(CDRH_distance_temp)
CDRH_distance[CDRH_distance > hamming_distance] <- NA  #set hamming distance limit
# Change axis names
colnames(CDRH_distance) <- matrix_sequences
CDRH_distance$cdrh3 <- colnames(CDRH_distance)
CDRH_distance <- CDRH_distance[,c(which(colnames(CDRH_distance)=="cdrh3"),which(colnames(CDRH_distance)!="cdrh3"))]
# Column name to first row and remove labels for duplicates
CDRH_distance <- rbind(names(CDRH_distance), CDRH_distance)
CDRH_distance[1,] <- gsub(x = CDRH_distance[1,], pattern = "\\..*", replacement = "")  
# Define clusters
CDRH_distance_cluster <- CDRH_distance
cluster_info <- as.data.frame(CDRH_distance_cluster$cdrh3)
cluster_info <- as.data.frame(cluster_info[-1,])
colnames(cluster_info)[which(names(cluster_info) == "cluster_info[-1, ]")] <- "cdrh3" # Change column name to be readable
cluster_info$cluster <- cluster_info$cdrh3
iteration_count <- 1 # Assign initial iteration number
hd <- hamming_distance + 1

repeat
{
  group_seed <- CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[,2] < hd)]
  for (s in group_seed)
  {
    # Designate layer object
    n <- 1
    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[s]] < hd)])
    # Assign cluster number and remove from dataframe
    cluster_info$cluster[cluster_info$cluster %in% s] <- iteration_count
    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == s]
    for (t in get(paste0("group_layer_",n)))
    {
      # Assign cluster number
      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
      # Assign next layer
      n <- n + 1
      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
      if (length(get(paste0("group_layer_",n))) == 0)
      {
        n <- n - 1
        next
      }
      # Remove assigned CDR from dataframe
      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
      for (t in get(paste0("group_layer_",n)))
      {
        # Assign cluster number
        cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
        # Assign next layer
        n <- n + 1
        assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
        if (length(get(paste0("group_layer_",n))) == 0)
        {
          n <- n - 1
          next
        }
        # Remove assigned CDR from dataframe
        CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
        for (t in get(paste0("group_layer_",n)))
        {
          # Assign cluster number
          cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
          # Assign next layer
          n <- n + 1
          assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
          if (length(get(paste0("group_layer_",n))) == 0)
          {
            n <- n - 1
            next
          }
          # Remove assigned CDR from dataframe
          CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
          for (t in get(paste0("group_layer_",n)))
          {
            # Assign cluster number
            cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
            # Assign next layer
            n <- n + 1
            assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
            if (length(get(paste0("group_layer_",n))) == 0)
            {
              n <- n - 1
              next
            }
            # Remove assigned CDR from dataframe
            CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
            for (t in get(paste0("group_layer_",n)))
            {
              # Assign cluster number
              cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
              # Assign next layer
              n <- n + 1
              assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
              if (length(get(paste0("group_layer_",n))) == 0)
              {
                n <- n - 1
                next
              }
              # Remove assigned CDR from dataframe
              CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
              for (t in get(paste0("group_layer_",n)))
              {
                # Assign cluster number
                cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                # Assign next layer
                n <- n + 1
                assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                if (length(get(paste0("group_layer_",n))) == 0)
                {
                  n <- n - 1
                  next
                }
                # Remove assigned CDR from dataframe
                CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                for (t in get(paste0("group_layer_",n)))
                {
                  # Assign cluster number
                  cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                  # Assign next layer
                  n <- n + 1
                  assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                  if (length(get(paste0("group_layer_",n))) == 0)
                  {
                    n <- n - 1
                    next
                  }
                  # Remove assigned CDR from dataframe
                  CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                  for (t in get(paste0("group_layer_",n)))
                  {
                    # Assign cluster number
                    cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                    # Assign next layer
                    n <- n + 1
                    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                    if (length(get(paste0("group_layer_",n))) == 0)
                    {
                      n <- n - 1
                      next
                    }
                    # Remove assigned CDR from dataframe
                    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                    for (t in get(paste0("group_layer_",n)))
                    {
                      # Assign cluster number
                      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                      # Assign next layer
                      n <- n + 1
                      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                      if (length(get(paste0("group_layer_",n))) == 0)
                      {
                        n <- n - 1
                        next
                      }
                      # Remove assigned CDR from dataframe
                      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                      for (t in get(paste0("group_layer_",n)))
                      {
                        # Assign cluster number
                        cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                        # Assign next layer
                        n <- n + 1
                        assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                        if (length(get(paste0("group_layer_",n))) == 0)
                        {
                          n <- n - 1
                          next
                        }
                        # Remove assigned CDR from dataframe
                        CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                        for (t in get(paste0("group_layer_",n)))
                        {
                          # Assign cluster number
                          cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                          # Assign next layer
                          n <- n + 1
                          assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                          if (length(get(paste0("group_layer_",n))) == 0)
                          {
                            n <- n - 1
                            next
                          }
                          # Remove assigned CDR from dataframe
                          CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                          for (t in get(paste0("group_layer_",n)))
                          {
                            # Assign cluster number
                            cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                            # Assign next layer
                            n <- n + 1
                            assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                            if (length(get(paste0("group_layer_",n))) == 0)
                            {
                              n <- n - 1
                              next
                            }
                            # Remove assigned CDR from dataframe
                            CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                            for (t in get(paste0("group_layer_",n)))
                            {
                              # Assign cluster number
                              cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                              # Assign next layer
                              n <- n + 1
                              assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                              if (length(get(paste0("group_layer_",n))) == 0)
                              {
                                n <- n - 1
                                next
                              }
                              # Remove assigned CDR from dataframe
                              CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                              for (t in get(paste0("group_layer_",n)))
                              {
                                # Assign cluster number
                                cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                # Assign next layer
                                n <- n + 1
                                assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                if (length(get(paste0("group_layer_",n))) == 0)
                                {
                                  n <- n - 1
                                  next
                                }
                                # Remove assigned CDR from dataframe
                                CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                for (t in get(paste0("group_layer_",n)))
                                {
                                  # Assign cluster number
                                  cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                  # Assign next layer
                                  n <- n + 1
                                  assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                  if (length(get(paste0("group_layer_",n))) == 0)
                                  {
                                    n <- n - 1
                                    next
                                  }
                                  # Remove assigned CDR from dataframe
                                  CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                  for (t in get(paste0("group_layer_",n)))
                                  {
                                    # Assign cluster number
                                    cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                    # Assign next layer
                                    n <- n + 1
                                    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                    if (length(get(paste0("group_layer_",n))) == 0)
                                    {
                                      n <- n - 1
                                      next
                                    }
                                    # Remove assigned CDR from dataframe
                                    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                    for (t in get(paste0("group_layer_",n)))
                                    {
                                      # Assign cluster number
                                      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                      # Assign next layer
                                      n <- n + 1
                                      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                      if (length(get(paste0("group_layer_",n))) == 0)
                                      {
                                        n <- n - 1
                                        next
                                      }
                                      # Remove assigned CDR from dataframe
                                      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                      for (t in get(paste0("group_layer_",n)))
                                      {
                                        # Assign cluster number
                                        cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                        # Assign next layer
                                        n <- n + 1
                                        assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                        if (length(get(paste0("group_layer_",n))) == 0)
                                        {
                                          n <- n - 1
                                          next
                                        }
                                        # Remove assigned CDR from dataframe
                                        CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                        for (t in get(paste0("group_layer_",n)))
                                        {
                                          # Assign cluster number
                                          cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                          # Assign next layer
                                          n <- n + 1
                                          assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                          if (length(get(paste0("group_layer_",n))) == 0)
                                          {
                                            n <- n - 1
                                            next
                                          }
                                          # Remove assigned CDR from dataframe
                                          CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                          for (t in get(paste0("group_layer_",n)))
                                          {
                                            # Assign cluster number
                                            cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                            # Assign next layer
                                            n <- n + 1
                                            assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                            if (length(get(paste0("group_layer_",n))) == 0)
                                            {
                                              n <- n - 1
                                              next
                                            }
                                            # Remove assigned CDR from dataframe
                                            CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                            for (t in get(paste0("group_layer_",n)))
                                            {
                                              # Assign cluster number
                                              cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                              # Assign next layer
                                              n <- n + 1
                                              assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                              if (length(get(paste0("group_layer_",n))) == 0)
                                              {
                                                n <- n - 1
                                                next
                                              }
                                              # Remove assigned CDR from dataframe
                                              CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                              for (t in get(paste0("group_layer_",n)))
                                              {
                                                # Assign cluster number
                                                cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                # Assign next layer
                                                n <- n + 1
                                                assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                if (length(get(paste0("group_layer_",n))) == 0)
                                                {
                                                  n <- n - 1
                                                  next
                                                }
                                                # Remove assigned CDR from dataframe
                                                CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                for (t in get(paste0("group_layer_",n)))
                                                {
                                                  # Assign cluster number
                                                  cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                  # Assign next layer
                                                  n <- n + 1
                                                  assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                  if (length(get(paste0("group_layer_",n))) == 0)
                                                  {
                                                    n <- n - 1
                                                    next
                                                  }
                                                  # Remove assigned CDR from dataframe
                                                  CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                  for (t in get(paste0("group_layer_",n)))
                                                  {
                                                    # Assign cluster number
                                                    cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                    # Assign next layer
                                                    n <- n + 1
                                                    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                    if (length(get(paste0("group_layer_",n))) == 0)
                                                    {
                                                      n <- n - 1
                                                      next
                                                    }
                                                    # Remove assigned CDR from dataframe
                                                    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                    for (t in get(paste0("group_layer_",n)))
                                                    {
                                                      # Assign cluster number
                                                      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                      # Assign next layer
                                                      n <- n + 1
                                                      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                      if (length(get(paste0("group_layer_",n))) == 0)
                                                      {
                                                        n <- n - 1
                                                        next
                                                      }
                                                      # Remove assigned CDR from dataframe
                                                      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                      for (t in get(paste0("group_layer_",n)))
                                                      {
                                                        # Assign cluster number
                                                        cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                        # Assign next layer
                                                        n <- n + 1
                                                        assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                        if (length(get(paste0("group_layer_",n))) == 0)
                                                        {
                                                          n <- n - 1
                                                          next
                                                        }
                                                        # Remove assigned CDR from dataframe
                                                        CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                        for (t in get(paste0("group_layer_",n)))
                                                        {
                                                          # Assign cluster number
                                                          cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                          # Assign next layer
                                                          n <- n + 1
                                                          assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                          if (length(get(paste0("group_layer_",n))) == 0)
                                                          {
                                                            n <- n - 1
                                                            next
                                                          }
                                                          # Remove assigned CDR from dataframe
                                                          CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                          for (t in get(paste0("group_layer_",n)))
                                                          {
                                                            # Assign cluster number
                                                            cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                            # Assign next layer
                                                            n <- n + 1
                                                            assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                            if (length(get(paste0("group_layer_",n))) == 0)
                                                            {
                                                              n <- n - 1
                                                              next
                                                            }
                                                            # Remove assigned CDR from dataframe
                                                            CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                            for (t in get(paste0("group_layer_",n)))
                                                            {
                                                              # Assign cluster number
                                                              cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                              # Assign next layer
                                                              n <- n + 1
                                                              assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                              if (length(get(paste0("group_layer_",n))) == 0)
                                                              {
                                                                n <- n - 1
                                                                next
                                                              }
                                                              # Remove assigned CDR from dataframe
                                                              CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                              for (t in get(paste0("group_layer_",n)))
                                                              {
                                                                # Assign cluster number
                                                                cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                # Assign next layer
                                                                n <- n + 1
                                                                assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                if (length(get(paste0("group_layer_",n))) == 0)
                                                                {
                                                                  n <- n - 1
                                                                  next
                                                                }
                                                                # Remove assigned CDR from dataframe
                                                                CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                for (t in get(paste0("group_layer_",n)))
                                                                {
                                                                  # Assign cluster number
                                                                  cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                  # Assign next layer
                                                                  n <- n + 1
                                                                  assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                  if (length(get(paste0("group_layer_",n))) == 0)
                                                                  {
                                                                    n <- n - 1
                                                                    next
                                                                  }
                                                                  # Remove assigned CDR from dataframe
                                                                  CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                  for (t in get(paste0("group_layer_",n)))
                                                                  {
                                                                    # Assign cluster number
                                                                    cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                    # Assign next layer
                                                                    n <- n + 1
                                                                    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                    if (length(get(paste0("group_layer_",n))) == 0)
                                                                    {
                                                                      n <- n - 1
                                                                      next
                                                                    }
                                                                    # Remove assigned CDR from dataframe
                                                                    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                    for (t in get(paste0("group_layer_",n)))
                                                                    {
                                                                      # Assign cluster number
                                                                      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                      # Assign next layer
                                                                      n <- n + 1
                                                                      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                      if (length(get(paste0("group_layer_",n))) == 0)
                                                                      {
                                                                        n <- n - 1
                                                                        next
                                                                      }
                                                                      # Remove assigned CDR from dataframe
                                                                      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                      for (t in get(paste0("group_layer_",n)))
                                                                      {
                                                                        # Assign cluster number
                                                                        cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                        # Assign next layer
                                                                        n <- n + 1
                                                                        assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                        if (length(get(paste0("group_layer_",n))) == 0)
                                                                        {
                                                                          n <- n - 1
                                                                          next
                                                                        }
                                                                        # Remove assigned CDR from dataframe
                                                                        CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                        for (t in get(paste0("group_layer_",n)))
                                                                        {
                                                                          # Assign cluster number
                                                                          cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                          # Assign next layer
                                                                          n <- n + 1
                                                                          assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                          if (length(get(paste0("group_layer_",n))) == 0)
                                                                          {
                                                                            n <- n - 1
                                                                            next
                                                                          }
                                                                          # Remove assigned CDR from dataframe
                                                                          CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                          for (t in get(paste0("group_layer_",n)))
                                                                          {
                                                                            # Assign cluster number
                                                                            cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                            # Assign next layer
                                                                            n <- n + 1
                                                                            assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                            if (length(get(paste0("group_layer_",n))) == 0)
                                                                            {
                                                                              n <- n - 1
                                                                              next
                                                                            }
                                                                            # Remove assigned CDR from dataframe
                                                                            CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                            for (t in get(paste0("group_layer_",n)))
                                                                            {
                                                                              # Assign cluster number
                                                                              cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                              # Assign next layer
                                                                              n <- n + 1
                                                                              assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                              if (length(get(paste0("group_layer_",n))) == 0)
                                                                              {
                                                                                n <- n - 1
                                                                                next
                                                                              }
                                                                              # Remove assigned CDR from dataframe
                                                                              CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                              for (t in get(paste0("group_layer_",n)))
                                                                              {
                                                                                # Assign cluster number
                                                                                cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                                # Assign next layer
                                                                                n <- n + 1
                                                                                assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                                if (length(get(paste0("group_layer_",n))) == 0)
                                                                                {
                                                                                  n <- n - 1
                                                                                  next
                                                                                }
                                                                                # Remove assigned CDR from dataframe
                                                                                CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                                for (t in get(paste0("group_layer_",n)))
                                                                                {
                                                                                  # Assign cluster number
                                                                                  cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                                  # Assign next layer
                                                                                  n <- n + 1
                                                                                  assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                                  if (length(get(paste0("group_layer_",n))) == 0)
                                                                                  {
                                                                                    n <- n - 1
                                                                                    next
                                                                                  }
                                                                                  # Remove assigned CDR from dataframe
                                                                                  CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                                  for (t in get(paste0("group_layer_",n)))
                                                                                  {
                                                                                    # Assign cluster number
                                                                                    cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                                    # Assign next layer
                                                                                    n <- n + 1
                                                                                    assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                                    if (length(get(paste0("group_layer_",n))) == 0)
                                                                                    {
                                                                                      n <- n - 1
                                                                                      next
                                                                                    }
                                                                                    # Remove assigned CDR from dataframe
                                                                                    CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                                    for (t in get(paste0("group_layer_",n)))
                                                                                    {
                                                                                      # Assign cluster number
                                                                                      cluster_info$cluster[cluster_info$cluster %in% t] <- iteration_count
                                                                                      # Assign next layer
                                                                                      n <- n + 1
                                                                                      assign(paste0("group_layer_",n), CDRH_distance_cluster$cdrh3[which(CDRH_distance_cluster[[t]] < hd)])
                                                                                      if (length(get(paste0("group_layer_",n))) == 0)
                                                                                      {
                                                                                        n <- n - 1
                                                                                        next
                                                                                      }
                                                                                      # Remove assigned CDR from dataframe
                                                                                      CDRH_distance_cluster <- CDRH_distance_cluster[,!CDRH_distance_cluster[1,] == t]
                                                                                    }
                                                                                  }
                                                                                }
                                                                              }
                                                                            }
                                                                          }
                                                                        }
                                                                      }
                                                                    }
                                                                  }
                                                                }
                                                              }
                                                            }
                                                          }
                                                        }
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  iteration_count <- iteration_count + 1
  if (ncol(CDRH_distance_cluster) == 1) {break}
}

# Save data
save(cluster_info, file = paste0("Directory/CDRH3_cluster_", IGHV, "_hamming_", hamming_distance,".Rdata"))





