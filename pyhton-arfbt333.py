import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch
import os
import glob
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


class CTMaterialAnalyzer:
    def __init__(self):
        """
        Initialize CT Material Analyzer
        """
        # Image dimensions
        self.width = 105
        self.height = 105

        # Material color mapping (consistent with C code)
        # 注意：将Unknown放在第一个位置，这样-1会映射到Unknown的灰色
        self.material_colors = {
            'Unknown': (200 / 255, 200 / 255, 200 / 255),  # Gray - 放在第一位
            'LiMetal': (255 / 255, 0 / 255, 0 / 255),  # Red
            'MnO2': (0 / 255, 0 / 255, 255 / 255),  # Blue
            'Li0.1MnO2': (0 / 255, 50 / 255, 255 / 255),  # Light blue
            'Li0.2MnO2': (0 / 255, 100 / 255, 255 / 255),  # Blue-green
            'Li0.3MnO2': (0 / 255, 150 / 255, 255 / 255),  # Cyan-blue
            'Li0.4MnO2': (0 / 255, 200 / 255, 255 / 255),  # Cyan-green
            'Li0.5MnO2': (0 / 255, 255 / 255, 255 / 255),  # Cyan
            'Li0.6MnO2': (0 / 255, 255 / 255, 200 / 255),  # Cyan-yellow
            'Li0.7MnO2': (0 / 255, 255 / 255, 150 / 255),  # Yellow-green
            'Li0.8MnO2': (0 / 255, 255 / 255, 100 / 255),  # Orange-yellow
            'Li0.9MnO2': (0 / 255, 255 / 255, 50 / 255),  # Orange
            'LiMnO2': (0 / 255, 255 / 255, 0 / 255),  # Green
        }

        # Define region coordinates (consistent with C code)
        self.regions = {
            'fresh': {
                'material': {'min_x': 44, 'max_x': 65, 'min_y': 22, 'max_y': 85},
                'cathode': {'min_x': 47, 'max_x': 61, 'min_y': 22, 'max_y': 85},
                'anode': {'min_x': 63, 'max_x': 65, 'min_y': 22, 'max_y': 85}
            },
            'discharged': {
                'material': {'min_x': 45, 'max_x': 67, 'min_y': 30, 'max_y': 83},
                'cathode': {'min_x': 49, 'max_x': 60, 'min_y': 31, 'max_y': 82},
                'anode': {'min_x': 63, 'max_x': 65, 'min_y': 22, 'max_y': 85}
            }
        }

        # Create colormaps
        self.alpha_cmap = LinearSegmentedColormap.from_list(
            'alpha_cmap', ['blue', 'cyan', 'green', 'yellow', 'red']
        )

        self.soc_cmap = LinearSegmentedColormap.from_list(
            'soc_cmap', ['darkblue', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red']
        )

        print("CT Material Analyzer initialized")

    def load_material_data(self, csv_file):
        """
        Load material identification CSV file

        Parameters:
        csv_file: str, CSV file path
        """
        try:
            # Try different encodings
            encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            for encoding in encodings:
                try:
                    self.df = pd.read_csv(csv_file, encoding=encoding)
                    print(f"Successfully loaded with {encoding} encoding")
                    break
                except UnicodeDecodeError:
                    continue
            else:
                # If all encodings fail, try without specifying encoding
                self.df = pd.read_csv(csv_file)

            # Determine state from filename
            filename_lower = csv_file.lower()
            if 'fresh' in filename_lower:
                self.state = 'fresh'
            elif 'discharged' in filename_lower:
                self.state = 'discharged'
            else:
                self.state = 'unknown'
                print("Warning: Cannot determine battery state from filename")
                self.state = 'fresh'

            print(f"Data loaded successfully: {csv_file}")
            print(f"Data shape: {self.df.shape}")
            print(f"Detected state: {self.state}")
            print(f"Data columns: {self.df.columns.tolist()}")

            # For discharged state, remove lithium metal identification
            if self.state == 'discharged' and 'Material' in self.df.columns:
                self.df.loc[self.df['Material'] == 'LiMetal', 'Material'] = 'Unknown'
                print("Discharged state: Lithium metal reclassified as Unknown")

            # For fresh state, only keep lithium metal in anode region
            elif self.state == 'fresh' and 'Material' in self.df.columns:
                anode_region = self.regions[self.state]['anode']

                # 创建掩码：识别LiMetal且在阳极区域外的点
                mask_li_metal = self.df['Material'] == 'LiMetal'
                mask_not_in_anode = ~(
                        (self.df['X'] >= anode_region['min_x']) &
                        (self.df['X'] <= anode_region['max_x']) &
                        (self.df['Y'] >= anode_region['min_y']) &
                        (self.df['Y'] <= anode_region['max_y'])
                )

                # 将阳极区域外的LiMetal重新分类为Unknown
                li_metal_outside_anode = mask_li_metal & mask_not_in_anode
                if li_metal_outside_anode.any():
                    count = li_metal_outside_anode.sum()
                    self.df.loc[li_metal_outside_anode, 'Material'] = 'Unknown'
                    print(f"Fresh state: {count} LiMetal pixels outside anode region reclassified as Unknown")
                else:
                    print("Fresh state: All LiMetal pixels are within anode region")

            # Check required columns
            required_cols = ['X', 'Y', 'Material', 'Pe', 'Zeff', 'SOC']
            missing_cols = [col for col in required_cols if col not in self.df.columns]
            if missing_cols:
                print(f"Warning: Missing columns {missing_cols}")

            return True

        except Exception as e:
            print(f"Failed to load data: {e}")
            return False

    def get_region_data(self, region_type='material'):
        """
        Get data for specified region

        Parameters:
        region_type: str, region type ('material', 'cathode', or 'anode')
        """
        if self.state not in self.regions:
            return self.df

        if region_type not in self.regions[self.state]:
            print(f"Warning: Region type '{region_type}' not defined for state '{self.state}'")
            return self.df

        region = self.regions[self.state][region_type]

        # Ensure columns exist
        if 'X' not in self.df.columns or 'Y' not in self.df.columns:
            print("Error: Data missing X or Y columns")
            return pd.DataFrame()

        mask = (
                (self.df['X'] >= region['min_x']) &
                (self.df['X'] <= region['max_x']) &
                (self.df['Y'] >= region['min_y']) &
                (self.df['Y'] <= region['max_y'])
        )

        return self.df[mask]

    def create_material_distribution_map(self, save_path=None):
        """
        Create material distribution map

        Parameters:
        save_path: str, save path
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        # Material distribution map - 初始化为Unknown（灰色）
        material_map = np.full((self.height, self.width), 0, dtype=int)  # 0对应Unknown
        material_names = list(self.material_colors.keys())

        # Check required columns
        required_cols = ['X', 'Y', 'Material']
        missing_cols = [col for col in required_cols if col not in self.df.columns]
        if missing_cols:
            print(f"Error: Missing required columns {missing_cols}")
            return

        # 获取材料区域的数据
        material_region_data = self.get_region_data('material')

        # 只绘制材料区域内的点
        for _, row in material_region_data.iterrows():
            try:
                x, y = int(row['X']), int(row['Y'])
                material = str(row['Material'])
                if material in material_names:
                    material_map[y, x] = material_names.index(material)
                else:
                    material_map[y, x] = 0  # Unknown
            except Exception:
                continue

        # Create colormap
        colors = list(self.material_colors.values())
        cmap = mcolors.ListedColormap(colors)
        bounds = range(len(material_names) + 1)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        im1 = ax1.imshow(material_map, cmap=cmap, norm=norm, origin='lower')
        ax1.set_title(f'Material Distribution - {self.state.capitalize()}', fontsize=14, fontweight='bold')
        ax1.set_xlabel('X (pixel)')
        ax1.set_ylabel('Y (pixel)')

        # Highlight cathode region - 保留矩形框，移除文字标注
        if self.state in self.regions:
            region = self.regions[self.state]['cathode']
            rect = Rectangle((region['min_x'], region['min_y']),
                             region['max_x'] - region['min_x'],
                             region['max_y'] - region['min_y'],
                             linewidth=2, edgecolor='red', facecolor='none')
            ax1.add_patch(rect)
            # 移除文字标注: ax1.text(region['min_x'], region['min_y'] - 5, 'Cathode Region',
            #             color='red', fontweight='bold')

            # Highlight anode region (only for fresh state) - 保留矩形框，移除文字标注
            if self.state == 'fresh':
                anode_region = self.regions[self.state]['anode']
                rect_anode = Rectangle((anode_region['min_x'], anode_region['min_y']),
                                       anode_region['max_x'] - anode_region['min_x'],
                                       anode_region['max_y'] - anode_region['min_y'],
                                       linewidth=2, edgecolor='orange', facecolor='none', linestyle='--')
                ax1.add_patch(rect_anode)
                # 移除文字标注: ax1.text(anode_region['min_x'], anode_region['min_y'] - 5, 'Anode Region',
                #          color='orange', fontweight='bold')

            # Highlight material region boundary - 保留矩形框，移除文字标注
            material_region = self.regions[self.state]['material']
            rect_material = Rectangle((material_region['min_x'], material_region['min_y']),
                                      material_region['max_x'] - material_region['min_x'],
                                      material_region['max_y'] - material_region['min_y'],
                                      linewidth=1, edgecolor='black', facecolor='none', linestyle=':')
            ax1.add_patch(rect_material)

        # 创建图例 - 跳过放电状态的LiMetal
        legend_elements = []
        for material, color in self.material_colors.items():
            if self.state == 'discharged' and material == 'LiMetal':
                continue
            legend_elements.append(Patch(facecolor=color, label=material))

        ax1.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

        # SOC distribution map - 也限制在材料区域内
        soc_map = np.zeros((self.height, self.width))
        if 'SOC' in material_region_data.columns:
            for _, row in material_region_data.iterrows():
                try:
                    x, y = int(row['X']), int(row['Y'])
                    soc_map[y, x] = float(row['SOC'])
                except Exception:
                    continue
        else:
            print("Warning: Data missing SOC column")

        im2 = ax2.imshow(soc_map, cmap=self.soc_cmap, origin='lower', vmin=0, vmax=100)
        ax2.set_title(f'SOC Distribution - {self.state.capitalize()}', fontsize=14, fontweight='bold')
        ax2.set_xlabel('X (pixel)')
        ax2.set_ylabel('Y (pixel)')

        # Highlight cathode region in SOC map - 保留矩形框，移除文字标注
        if self.state in self.regions:
            region = self.regions[self.state]['cathode']
            rect = Rectangle((region['min_x'], region['min_y']),
                             region['max_x'] - region['min_x'],
                             region['max_y'] - region['min_y'],
                             linewidth=2, edgecolor='red', facecolor='none')
            ax2.add_patch(rect)

            # Highlight anode region (only for fresh state) - 保留矩形框，移除文字标注
            if self.state == 'fresh':
                anode_region = self.regions[self.state]['anode']
                rect_anode = Rectangle((anode_region['min_x'], anode_region['min_y']),
                                       anode_region['max_x'] - anode_region['min_x'],
                                       anode_region['max_y'] - anode_region['min_y'],
                                       linewidth=2, edgecolor='orange', facecolor='none', linestyle='--')
                ax2.add_patch(rect_anode)

            # Highlight material region boundary - 保留矩形框，移除文字标注
            material_region = self.regions[self.state]['material']
            rect_material = Rectangle((material_region['min_x'], material_region['min_y']),
                                      material_region['max_x'] - material_region['min_x'],
                                      material_region['max_y'] - material_region['min_y'],
                                      linewidth=1, edgecolor='black', facecolor='none', linestyle=':')
            ax2.add_patch(rect_material)

        plt.colorbar(im2, ax=ax2, label='SOC (%)')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Material distribution map saved: {save_path}")

        plt.show()

    def create_pe_zeff_analysis(self, save_path=None):
        """
        Create Pe-Zeff analysis plot

        Parameters:
        save_path: str, save path
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        if isinstance(axes, np.ndarray):
            ax1, ax2, ax3 = axes
        else:
            ax1 = axes
            ax2, ax3 = None, None

        # Check required columns
        required_cols = ['Pe', 'Zeff', 'Material']
        missing_cols = [col for col in required_cols if col not in self.df.columns]
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols}")
            plt.close(fig)
            return

        # 只使用材料区域内的数据进行绘制
        material_region_data = self.get_region_data('material')

        # Scatter plot colored by material type
        for material, color in self.material_colors.items():
            if material == 'Unknown':
                continue
            # Skip LiMetal for discharged state
            if self.state == 'discharged' and material == 'LiMetal':
                continue

            mask = material_region_data['Material'] == material
            if mask.any():
                ax1.scatter(material_region_data.loc[mask, 'Pe'], material_region_data.loc[mask, 'Zeff'],
                            c=[color], label=material, alpha=0.7, s=10)

        ax1.set_xlabel('Pe (e/cm³)')
        ax1.set_ylabel('Zeff')
        ax1.set_title(f'Pe vs Zeff (by Material) - {self.state.capitalize()}')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')

        # Scatter plot colored by SOC
        if 'SOC' in material_region_data.columns:
            scatter = ax2.scatter(material_region_data['Pe'], material_region_data['Zeff'],
                                  c=material_region_data['SOC'], cmap=self.soc_cmap,
                                  alpha=0.7, s=10, vmin=0, vmax=100)
            ax2.set_xlabel('Pe (e/cm³)')
            ax2.set_ylabel('Zeff')
            ax2.set_title(f'Pe vs Zeff (by SOC) - {self.state.capitalize()}')
            plt.colorbar(scatter, ax=ax2, label='SOC (%)')
            ax2.grid(True, alpha=0.3)
            ax2.set_xscale('log')
        else:
            ax2.set_visible(False)

        # Cathode region analysis (discharged state)
        if self.state == 'discharged' and 'Alpha' in self.df.columns and ax3 is not None:
            cathode_df = self.get_region_data('cathode')
            if len(cathode_df) > 0 and 'Alpha' in cathode_df.columns:
                valid_cathode = cathode_df[cathode_df['Alpha'] >= 0]
                if len(valid_cathode) > 0:
                    sc = ax3.scatter(valid_cathode['Pe'], valid_cathode['Zeff'],
                                     c=valid_cathode['Alpha'], cmap=self.alpha_cmap,
                                     alpha=0.7, s=20, vmin=0, vmax=1)
                    ax3.set_xlabel('Pe (e/cm³)')
                    ax3.set_ylabel('Zeff')
                    ax3.set_title('Pe vs Zeff in Cathode (by α)')
                    ax3.grid(True, alpha=0.3)
                    ax3.set_xscale('log')
                    plt.colorbar(sc, ax=ax3, label='α value')
                else:
                    ax3.set_visible(False)
            else:
                ax3.set_visible(False)
        elif ax3 is not None:
            ax3.set_visible(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Pe-Zeff analysis plot saved: {save_path}")

        plt.show()

    def create_alpha_beta_analysis(self, save_path=None):
        """
        Create α-β analysis plot (discharged state only)

        Parameters:
        save_path: str, save path
        """
        if self.state != 'discharged':
            print("α-β analysis only for discharged state")
            return

        if 'Alpha' not in self.df.columns or 'Beta' not in self.df.columns:
            print("Data missing Alpha or Beta columns")
            return

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

        # 只使用材料区域内的数据进行绘制
        material_region_data = self.get_region_data('material')

        # α value distribution map - 只绘制材料区域内的点
        alpha_map = np.full((self.height, self.width), np.nan, dtype=float)  # 使用NaN表示无数据
        for _, row in material_region_data.iterrows():
            try:
                x, y = int(row['X']), int(row['Y'])
                alpha = float(row['Alpha'])
                if alpha >= 0:
                    alpha_map[y, x] = alpha
            except Exception:
                continue

        im1 = ax1.imshow(alpha_map, cmap=self.alpha_cmap, origin='lower', vmin=0, vmax=1)
        ax1.set_title('α(Li Content) Distribution', fontsize=14, fontweight='bold')
        ax1.set_xlabel('X (pixel)')
        ax1.set_ylabel('Y (pixel)')

        # Highlight cathode region - 保留矩形框，移除文字标注
        if self.state in self.regions:
            region = self.regions[self.state]['cathode']
            rect = Rectangle((region['min_x'], region['min_y']),
                             region['max_x'] - region['min_x'],
                             region['max_y'] - region['min_y'],
                             linewidth=2, edgecolor='red', facecolor='none')
            ax1.add_patch(rect)

            # Highlight material region boundary - 保留矩形框，移除文字标注
            material_region = self.regions[self.state]['material']
            rect_material = Rectangle((material_region['min_x'], material_region['min_y']),
                                      material_region['max_x'] - material_region['min_x'],
                                      material_region['max_y'] - material_region['min_y'],
                                      linewidth=1, edgecolor='black', facecolor='none', linestyle=':')
            ax1.add_patch(rect_material)

        plt.colorbar(im1, ax=ax1, label='α value')

        # β value distribution map - 只绘制材料区域内的点
        beta_map = np.full((self.height, self.width), np.nan, dtype=float)  # 使用NaN表示无数据
        for _, row in material_region_data.iterrows():
            try:
                x, y = int(row['X']), int(row['Y'])
                beta = float(row['Beta'])
                if beta >= 0:
                    beta_map[y, x] = beta
            except Exception:
                continue

        im2 = ax2.imshow(beta_map, cmap=self.alpha_cmap, origin='lower', vmin=0, vmax=1)
        ax2.set_title('β(MnO₂ Content) Distribution', fontsize=14, fontweight='bold')
        ax2.set_xlabel('X (pixel)')
        ax2.set_ylabel('Y (pixel)')

        # Highlight cathode region - 保留矩形框，移除文字标注
        if self.state in self.regions:
            rect = Rectangle((region['min_x'], region['min_y']),
                             region['max_x'] - region['min_x'],
                             region['max_y'] - region['min_y'],
                             linewidth=2, edgecolor='red', facecolor='none')
            ax2.add_patch(rect)

            # Highlight material region boundary - 保留矩形框，移除文字标注
            rect_material = Rectangle((material_region['min_x'], material_region['min_y']),
                                      material_region['max_x'] - material_region['min_x'],
                                      material_region['max_y'] - material_region['min_y'],
                                      linewidth=1, edgecolor='black', facecolor='none', linestyle=':')
            ax2.add_patch(rect_material)

        plt.colorbar(im2, ax=ax2, label='β value')

        # α-β scatter plot
        cathode_df = self.get_region_data('cathode')
        if len(cathode_df) > 0:
            if 'Alpha' in cathode_df.columns and 'Beta' in cathode_df.columns:
                valid_cathode = cathode_df[(cathode_df['Alpha'] >= 0) & (cathode_df['Beta'] >= 0)]
                if len(valid_cathode) > 0 and 'SOC' in valid_cathode.columns:
                    sc = ax3.scatter(valid_cathode['Alpha'], valid_cathode['Beta'],
                                     c=valid_cathode['SOC'], cmap=self.soc_cmap,
                                     alpha=0.7, s=20, vmin=0, vmax=100)
                    ax3.set_xlabel('α value (Li content)')
                    ax3.set_ylabel('β value (MnO₂ content)')
                    ax3.set_title('α vs β in Cathode Region (colored by SOC)')
                    ax3.grid(True, alpha=0.3)

                    # Add diagonal line (α + β = 1)
                    ax3.plot([0, 1], [1, 0], 'r--', alpha=0.5, label='α + β = 1')
                    ax3.legend()

                    plt.colorbar(sc, ax=ax3, label='SOC (%)')
                else:
                    ax3.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax3.transAxes)
            else:
                ax3.text(0.5, 0.5, 'Missing Alpha or Beta columns', ha='center', va='center', transform=ax3.transAxes)
        else:
            ax3.text(0.5, 0.5, 'No cathode region data', ha='center', va='center', transform=ax3.transAxes)

        # α vs SOC relationship
        if len(cathode_df) > 0 and 'Alpha' in cathode_df.columns and 'SOC' in cathode_df.columns:
            valid_cathode = cathode_df[cathode_df['Alpha'] >= 0]
            if len(valid_cathode) > 0:
                sc = ax4.scatter(valid_cathode['Alpha'], valid_cathode['SOC'],
                                 c=valid_cathode['Alpha'], cmap=self.alpha_cmap,
                                 alpha=0.7, s=20, vmin=0, vmax=1)
                ax4.set_xlabel('α value (Li content)')
                ax4.set_ylabel('SOC (%)')
                ax4.set_title('α vs SOC in Cathode Region')
                ax4.grid(True, alpha=0.3)

                # Add trend line
                try:
                    z = np.polyfit(valid_cathode['Alpha'], valid_cathode['SOC'], 1)
                    p = np.poly1d(z)
                    ax4.plot(valid_cathode['Alpha'], p(valid_cathode['Alpha']), "r--", alpha=0.8)

                    correlation = np.corrcoef(valid_cathode['Alpha'], valid_cathode['SOC'])[0, 1]
                    ax4.text(0.05, 0.95, f'Correlation: {correlation:.3f}',
                             transform=ax4.transAxes, fontsize=12,
                             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))
                except Exception:
                    pass
            else:
                ax4.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax4.transAxes)
        else:
            ax4.text(0.5, 0.5, 'Missing Alpha or SOC columns', ha='center', va='center', transform=ax4.transAxes)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"α-β analysis plot saved: {save_path}")

        plt.show()

    def create_discharged_cathode_analysis(self, save_path=None):
        """
        Create detailed cathode region analysis for discharged state

        Parameters:
        save_path: str, save path
        """
        if self.state != 'discharged':
            print("Cathode analysis only for discharged state")
            return

        cathode_df = self.get_region_data('cathode')
        if len(cathode_df) == 0:
            print("No cathode region data")
            return

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

        # 1. Cathode region SOC distribution histogram
        if 'SOC' in cathode_df.columns:
            soc_bins = [0, 20, 40, 60, 80, 100]
            soc_labels = ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
            soc_dist = pd.cut(cathode_df['SOC'], bins=soc_bins, labels=soc_labels).value_counts().sort_index()

            colors = plt.cm.viridis(np.linspace(0, 1, len(soc_dist)))
            ax1.bar(soc_labels, soc_dist.values, color=colors)
            ax1.set_title('Cathode Region SOC Distribution')
            ax1.set_xlabel('SOC Range')
            ax1.set_ylabel('Pixel Count')
            ax1.tick_params(axis='x', rotation=45)

            for i, v in enumerate(soc_dist.values):
                ax1.text(i, v + max(soc_dist.values) * 0.01, str(v),
                         ha='center', va='bottom', fontweight='bold')
        else:
            ax1.text(0.5, 0.5, 'Missing SOC data', ha='center', va='center', transform=ax1.transAxes)

        # 2. Material composition by SOC range
        if 'SOC' in cathode_df.columns and 'Material' in cathode_df.columns:
            soc_ranges = {
                'Low SOC (0-33%)': (0, 33),
                'Mid SOC (33-66%)': (33, 66),
                'High SOC (66-100%)': (66, 100)
            }

            soc_range_data = []
            for range_name, (soc_min, soc_max) in soc_ranges.items():
                range_mask = (cathode_df['SOC'] >= soc_min) & (cathode_df['SOC'] <= soc_max)
                range_df = cathode_df[range_mask]
                if len(range_df) > 0:
                    material_counts = range_df['Material'].value_counts()
                    for material, count in material_counts.items():
                        soc_range_data.append({
                            'SOC Range': range_name,
                            'Material': material,
                            'Count': count,
                            'Percentage': (count / len(range_df)) * 100
                        })

            if soc_range_data:
                soc_range_df = pd.DataFrame(soc_range_data)
                pivot_df = soc_range_df.pivot(index='Material', columns='SOC Range', values='Percentage').fillna(0)
                pivot_df.plot(kind='bar', ax=ax2, color=['#ff9999', '#66b3ff', '#99ff99'])
                ax2.set_title('Material Composition by SOC Range')
                ax2.set_xlabel('Material')
                ax2.set_ylabel('Percentage (%)')
                ax2.legend(title='SOC Range')
                ax2.tick_params(axis='x', rotation=45)
            else:
                ax2.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax2.transAxes)
        else:
            ax2.text(0.5, 0.5, 'Missing SOC or Material data', ha='center', va='center', transform=ax2.transAxes)

        # 3. Lithiation degree distribution
        if 'Alpha' in cathode_df.columns:
            valid_alpha = cathode_df[cathode_df['Alpha'] >= 0]['Alpha']
            if len(valid_alpha) > 0:
                alpha_bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
                alpha_labels = ['0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5',
                                '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
                alpha_dist = pd.cut(valid_alpha, bins=alpha_bins, labels=alpha_labels).value_counts().sort_index()

                colors = plt.cm.plasma(np.linspace(0, 1, len(alpha_dist)))
                ax3.bar(alpha_labels, alpha_dist.values, color=colors)
                ax3.set_title('Lithiation Degree (α) Distribution in Cathode')
                ax3.set_xlabel('α Value Range')
                ax3.set_ylabel('Pixel Count')
                ax3.tick_params(axis='x', rotation=45)

                for i, v in enumerate(alpha_dist.values):
                    ax3.text(i, v + max(alpha_dist.values) * 0.01, str(v),
                             ha='center', va='bottom', fontweight='bold')
            else:
                ax3.text(0.5, 0.5, 'No valid α data', ha='center', va='center', transform=ax3.transAxes)
        else:
            ax3.text(0.5, 0.5, 'Missing Alpha data', ha='center', va='center', transform=ax3.transAxes)

        # 4. Cathode region SOC spatial distribution
        if 'SOC' in cathode_df.columns:
            # 只绘制阴极区域内的SOC分布
            cathode_soc_map = np.full((self.height, self.width), np.nan)
            for _, row in cathode_df.iterrows():
                try:
                    x, y = int(row['X']), int(row['Y'])
                    cathode_soc_map[y, x] = float(row['SOC'])
                except Exception:
                    continue

            im = ax4.imshow(cathode_soc_map, cmap=self.soc_cmap, origin='lower', vmin=0, vmax=100)
            ax4.set_title('Cathode Region SOC Spatial Distribution')
            ax4.set_xlabel('X (pixel)')
            ax4.set_ylabel('Y (pixel)')

            # 在图上标出阴极区域边界 - 保留矩形框，移除文字标注
            if self.state in self.regions:
                region = self.regions[self.state]['cathode']
                rect = Rectangle((region['min_x'], region['min_y']),
                                 region['max_x'] - region['min_x'],
                                 region['max_y'] - region['min_y'],
                                 linewidth=2, edgecolor='red', facecolor='none')
                ax4.add_patch(rect)

            plt.colorbar(im, ax=ax4, label='SOC (%)')
        else:
            ax4.text(0.5, 0.5, 'Missing SOC data', ha='center', va='center', transform=ax4.transAxes)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Discharged cathode analysis plot saved: {save_path}")

        plt.show()

    def generate_statistical_report(self):
        """
        Generate statistical analysis report
        """
        print("\n" + "=" * 70)
        print(f"CT Material Analysis Report - {self.state.upper()}")
        print("=" * 70)

        # 只分析材料区域内的数据
        material_region_data = self.get_region_data('material')

        # Basic statistics
        total_pixels = len(material_region_data)
        print(f"\nGeneral Statistics (Material Region):")
        print("-" * 40)
        print(f"  Total Pixels: {total_pixels}")

        if 'SOC' in material_region_data.columns:
            print(f"  SOC Range: {material_region_data['SOC'].min():.1f}% - {material_region_data['SOC'].max():.1f}%")
            print(f"  Average SOC: {material_region_data['SOC'].mean():.1f}%")
            print(f"  SOC Std Dev: {material_region_data['SOC'].std():.1f}%")

        # Material distribution statistics
        if 'Material' in material_region_data.columns:
            material_stats = material_region_data['Material'].value_counts()
            print(f"\nMaterial Distribution (Material Region):")
            print("-" * 40)
            for material, count in material_stats.items():
                percentage = (count / total_pixels) * 100
                print(f"  {material:<15}: {percentage:6.1f}% ({count:6d} pixels)")
        else:
            print(f"\nMaterial Distribution (Material Region):")
            print("-" * 40)
            print("  Warning: No Material column in data")

        # Cathode region analysis
        cathode_df = self.get_region_data('cathode')
        if len(cathode_df) > 0:
            cathode_total = len(cathode_df)
            print(f"\nCathode Region Analysis ({cathode_total} pixels):")
            print("-" * 40)

            if 'Material' in cathode_df.columns:
                cathode_stats = cathode_df['Material'].value_counts()
                for material, count in cathode_stats.items():
                    percentage = (count / cathode_total) * 100
                    print(f"  {material:<15}: {percentage:6.1f}% ({count:6d} pixels)")
            else:
                print("  Warning: No Material column in cathode data")

            if 'SOC' in cathode_df.columns:
                print(f"\nCathode Region SOC Statistics:")
                print("-" * 40)
                print(f"  SOC Range: {cathode_df['SOC'].min():.1f}% - {cathode_df['SOC'].max():.1f}%")
                print(f"  Average SOC: {cathode_df['SOC'].mean():.1f}%")
                print(f"  SOC Std Dev: {cathode_df['SOC'].std():.1f}%")
            else:
                print(f"\nCathode Region SOC Statistics:")
                print("-" * 40)
                print("  Warning: No SOC column in cathode data")

            # Lithiation degree analysis for discharged state
            if self.state == 'discharged' and 'Alpha' in cathode_df.columns:
                valid_alpha = cathode_df[cathode_df['Alpha'] >= 0]['Alpha']
                if len(valid_alpha) > 0:
                    avg_alpha = valid_alpha.mean()
                    std_alpha = valid_alpha.std()

                    print(f"\nLithiation Degree Analysis (Cathode Region):")
                    print("-" * 40)
                    print(f"  Average α: {avg_alpha:.4f} ± {std_alpha:.4f}")
                    print(f"  Average β: {1 - avg_alpha:.4f}")
                    print(f"  Approximate Formula: Li_{{{avg_alpha:.3f}}}MnO_{{{2 * (1 - avg_alpha):.3f}}}")

                    # Lithiation degree distribution
                    alpha_bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
                    alpha_labels = ['0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5',
                                    '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
                    alpha_dist = pd.cut(valid_alpha, bins=alpha_bins, labels=alpha_labels).value_counts().sort_index()

                    print(f"\n  Lithiation Degree Distribution:")
                    for bin_label, count in alpha_dist.items():
                        percentage = (count / len(valid_alpha)) * 100
                        print(f"    α = {bin_label}: {percentage:5.1f}% ({count:4d} pixels)")
                else:
                    print(f"\nLithiation Degree Analysis (Cathode Region):")
                    print("-" * 40)
                    print("  Warning: No valid α values")
            elif self.state == 'discharged':
                print(f"\nLithiation Degree Analysis (Cathode Region):")
                print("-" * 40)
                print("  Warning: No Alpha column in data")

        # Anode region analysis (only for fresh state)
        if self.state == 'fresh':
            anode_df = self.get_region_data('anode')
            if len(anode_df) > 0:
                anode_total = len(anode_df)
                print(f"\nAnode Region Analysis ({anode_total} pixels):")
                print("-" * 40)

                if 'Material' in anode_df.columns:
                    anode_stats = anode_df['Material'].value_counts()
                    for material, count in anode_stats.items():
                        percentage = (count / anode_total) * 100
                        print(f"  {material:<15}: {percentage:6.1f}% ({count:6d} pixels)")
                else:
                    print("  Warning: No Material column in anode data")

    def comprehensive_analysis(self, csv_file, output_dir='analysis_results'):
        """
        Perform comprehensive analysis

        Parameters:
        csv_file: str, CSV file path
        output_dir: str, output directory
        """
        if not self.load_material_data(csv_file):
            return

        # Create output directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        base_name = os.path.splitext(os.path.basename(csv_file))[0]

        print(f"\nStarting comprehensive analysis: {csv_file}")

        # 1. Material distribution map
        material_plot = os.path.join(output_dir, f'{base_name}_material_distribution.png')
        self.create_material_distribution_map(material_plot)

        # 2. Pe-Zeff analysis plot
        pe_zeff_plot = os.path.join(output_dir, f'{base_name}_pe_zeff_analysis.png')
        self.create_pe_zeff_analysis(pe_zeff_plot)

        # 3. Special analysis for discharged state
        if self.state == 'discharged':
            # α-β analysis plot
            alpha_beta_plot = os.path.join(output_dir, f'{base_name}_alpha_beta_analysis.png')
            self.create_alpha_beta_analysis(alpha_beta_plot)

            # Cathode region detailed analysis
            cathode_plot = os.path.join(output_dir, f'{base_name}_cathode_analysis.png')
            self.create_discharged_cathode_analysis(cathode_plot)

        # 4. Statistical report
        stats_file = os.path.join(output_dir, f'{base_name}_statistics.txt')
        with open(stats_file, 'w', encoding='utf-8') as f:
            import sys
            from io import StringIO

            old_stdout = sys.stdout
            new_stdout = StringIO()
            sys.stdout = new_stdout

            self.generate_statistical_report()

            output = new_stdout.getvalue()
            sys.stdout = old_stdout

            f.write(output)
            print(f"Statistical report saved: {stats_file}")

        # 5. Console output of statistical report
        self.generate_statistical_report()

        print(f"\nComprehensive analysis completed! Results saved in: {output_dir}")


def batch_analyze_all_files(result_dir='result', output_base='analysis_results'):
    """
    Batch analyze all files

    Parameters:
    result_dir: str, result files directory
    output_base: str, output base directory
    """
    if not os.path.exists(result_dir):
        print(f"Result directory does not exist: {result_dir}")
        return

    # Find all material distribution CSV files
    csv_pattern = os.path.join(result_dir, 'MaterialDistribution_*.csv')
    csv_files = glob.glob(csv_pattern)

    # If not found, try other patterns
    if not csv_files:
        csv_pattern = os.path.join(result_dir, '*.csv')
        csv_files = glob.glob(csv_pattern)

    if not csv_files:
        print(f"No CSV files found in {result_dir}")
        return

    print(f"Found {len(csv_files)} CSV files:")
    for file in csv_files:
        print(f"  - {file}")

    analyzer = CTMaterialAnalyzer()

    for csv_file in csv_files:
        # Determine state from filename
        filename_lower = csv_file.lower()
        if 'fresh' in filename_lower:
            state = 'fresh'
        elif 'discharged' in filename_lower:
            state = 'discharged'
        else:
            state = 'unknown'

        output_dir = os.path.join(output_base, state)

        print(f"\n{'=' * 60}")
        print(f"Analyzing file: {csv_file}")
        print(f"{'=' * 60}")

        analyzer.comprehensive_analysis(csv_file, output_dir)


# Usage example
if __name__ == "__main__":
    # Method 1: Automatic batch analysis of all files
    print("=== Automatic Batch Analysis ===")
    batch_analyze_all_files('result', 'analysis_results')

    # Method 2: Analyze specific files individually
    print("\n=== Individual File Analysis ===")
    analyzer = CTMaterialAnalyzer()

    # Analyze fresh state
    fresh_file = "result/MaterialDistribution_fresh.csv"
    if os.path.exists(fresh_file):
        analyzer.comprehensive_analysis(fresh_file, "analysis_results/fresh")
    else:
        print(f"File not found: {fresh_file}")

    # Analyze discharged state
    discharged_file = "result/MaterialDistribution_discharged.csv"
    if os.path.exists(discharged_file):
        analyzer.comprehensive_analysis(discharged_file, "analysis_results/discharged")
    else:
        print(f"File not found: {discharged_file}")
