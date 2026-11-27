// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

const GA_ID = 'G-2TSEMF47BK';
// Analytics only in production builds, so local previews don't pollute the data.
const analytics = process.env.NODE_ENV === 'production'
	? [
			{ tag: 'script', attrs: { async: true, src: `https://www.googletagmanager.com/gtag/js?id=${GA_ID}` } },
			{ tag: 'script', content: `window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag('js', new Date());gtag('config', '${GA_ID}');` },
		]
	: [];

export default defineConfig({
	site: 'https://odelaneau.github.io',
	base: '/GLIMPSE',
	redirects: {
		'/glimpse1': '/GLIMPSE/glimpse1/index.html',
		// URLs of the previous (Jekyll) site
		'/docs/installation': '/GLIMPSE/install/',
		'/docs/installation/static_binaries': '/GLIMPSE/install/',
		'/docs/installation/docker': '/GLIMPSE/install/docker/',
		'/docs/installation/build_from_source': '/GLIMPSE/install/build/',
		'/docs/installation/build_from_source/system_requirements': '/GLIMPSE/install/build/',
		'/docs/installation/build_from_source/required_libraries': '/GLIMPSE/install/build/',
		'/docs/installation/build_from_source/compile_glimpse2': '/GLIMPSE/install/build/',
		'/docs/documentation': '/GLIMPSE/tools/',
		'/docs/tutorials': '/GLIMPSE/tutorial/getting-started/',
		'/docs/tutorials/getting_started': '/GLIMPSE/tutorial/getting-started/',
		'/docs/tutorials/ukb_phasing': '/GLIMPSE/',
		'/docs/tutorials/ukb_wgs': '/GLIMPSE/',
		'/docs/quickstart': '/GLIMPSE/',
		'/docs/glossary': '/GLIMPSE/glossary/',
		'/CHANGELOG': '/GLIMPSE/changelog/',
		'/docs/documentation/chunk': '/GLIMPSE/tools/chunk/',
		'/docs/documentation/split_reference': '/GLIMPSE/tools/split_reference/',
		'/docs/documentation/phase': '/GLIMPSE/tools/phase/',
		'/docs/documentation/ligate': '/GLIMPSE/tools/ligate/',
		'/docs/documentation/concordance': '/GLIMPSE/tools/concordance/',
		'/docs/documentation/inspect': '/GLIMPSE/tools/inspect/',
	},
	integrations: [
		starlight({
			title: 'GLIMPSE2',
			description: 'Phasing and imputation of low-coverage whole-genome sequencing data.',
			expressiveCode: {
				themes: ['github-dark', 'github-light'],
			},
			components: { SocialIcons: './src/components/SocialIcons.astro' },
			favicon: '/favicon.ico',
			customCss: ['./src/styles/custom.css'],
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/odelaneau/GLIMPSE' },
			],
			editLink: {
				baseUrl: 'https://github.com/odelaneau/GLIMPSE/edit/master/docs/',
			},
			lastUpdated: true,
			head: analytics,
			sidebar: [
				{ label: 'Install', items: [
					{ label: 'Download and install', slug: 'install' },
					{ slug: 'install/docker' },
					{ slug: 'install/build' },
					{ slug: 'install/dependencies' },
				]},
				{ label: 'Tutorial', items: [{ slug: 'tutorial/getting-started' }] },
				{ label: 'Documentation', items: [
					{ label: 'Overview', slug: 'tools' },
					{ slug: 'tools/chunk' },
					{ slug: 'tools/split_reference' },
					{ slug: 'tools/phase' },
					{ slug: 'tools/ligate' },
					{ slug: 'tools/concordance' },
					{ slug: 'tools/inspect' },
				]},
				{ slug: 'glossary' },
				{ slug: 'changelog' },
				{ label: 'GLIMPSE1 (archive)', link: '/glimpse1' },
			],
		}),
	],
});
